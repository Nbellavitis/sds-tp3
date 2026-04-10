package ar.edu.itba.sds.tp3;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.FileTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

/**
 * Build compact analysis caches for System 1 simulation outputs.
 *
 * The cache stores the observables needed by the Python plotting scripts so
 * they do not have to parse the large raw simulation files repeatedly.
 */
public class SimulationAnalysisCacheBuilder {

    private static final double DEFAULT_DS = 0.2;
    private static final String CACHE_SUFFIX = ".analysis.v1.json.gz";

    private static final class SnapshotFrame {
        private final double time;
        private final double[] x;
        private final double[] y;
        private final double[] vx;
        private final double[] vy;
        private final boolean[] used;
        private final int usedCount;

        private SnapshotFrame(double time, double[] x, double[] y, double[] vx, double[] vy,
                              boolean[] used, int usedCount) {
            this.time = time;
            this.x = x;
            this.y = y;
            this.vx = vx;
            this.vy = vy;
            this.used = used;
            this.usedCount = usedCount;
        }
    }

    private static final class ShellGeometry {
        private final double[] edges;
        private final double[] centers;
        private final double[] areas;

        private ShellGeometry(double[] edges, double[] centers, double[] areas) {
            this.edges = edges;
            this.centers = centers;
            this.areas = areas;
        }
    }

    private static final class RunningProfiles {
        private final double[] countTime;
        private final double[] velocityTimeSum;
        private final double[] particleTime;
        private double totalTime;

        private RunningProfiles(int nShells) {
            this.countTime = new double[nShells];
            this.velocityTimeSum = new double[nShells];
            this.particleTime = new double[nShells];
            this.totalTime = 0.0;
        }
    }

    private static final class LinearFit {
        private final double slope;
        private final double intercept;
        private final double rSquared;

        private LinearFit(double slope, double intercept, double rSquared) {
            this.slope = slope;
            this.intercept = intercept;
            this.rSquared = rSquared;
        }
    }

    private static final class AnalysisResult {
        private final Path sourcePath;
        private final long sourceSize;
        private final FileTime sourceMtime;
        private final Map<String, Double> metadata;
        private final double[] cfcTimes;
        private final double[] cfcValues;
        private final LinearFit cfcFit;
        private final double[] fuTimes;
        private final double[] fuValues;
        private final double[] shellCenters;
        private final double[] rhoValues;
        private final double[] vValues;
        private final double[] jValues;

        private AnalysisResult(Path sourcePath, long sourceSize, FileTime sourceMtime,
                               Map<String, Double> metadata,
                               double[] cfcTimes, double[] cfcValues, LinearFit cfcFit,
                               double[] fuTimes, double[] fuValues,
                               double[] shellCenters, double[] rhoValues,
                               double[] vValues, double[] jValues) {
            this.sourcePath = sourcePath;
            this.sourceSize = sourceSize;
            this.sourceMtime = sourceMtime;
            this.metadata = metadata;
            this.cfcTimes = cfcTimes;
            this.cfcValues = cfcValues;
            this.cfcFit = cfcFit;
            this.fuTimes = fuTimes;
            this.fuValues = fuValues;
            this.shellCenters = shellCenters;
            this.rhoValues = rhoValues;
            this.vValues = vValues;
            this.jValues = jValues;
        }
    }

    public static void main(String[] args) throws IOException {
        Locale.setDefault(Locale.US);

        Path inputPath = null;
        Path cacheDir = null;
        double dS = DEFAULT_DS;

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--cache-dir" -> cacheDir = Paths.get(args[++i]);
                case "--dS" -> dS = Double.parseDouble(args[++i]);
                default -> {
                    if (inputPath == null) {
                        inputPath = Paths.get(args[i]);
                    } else {
                        throw new IllegalArgumentException("Unexpected argument: " + args[i]);
                    }
                }
            }
        }

        if (inputPath == null) {
            throw new IllegalArgumentException("Usage: SimulationAnalysisCacheBuilder <sim_file_or_directory> [--cache-dir <dir>] [--dS <value>]");
        }

        Path normalizedInput = inputPath.toAbsolutePath().normalize();
        Path resolvedCacheDir = (cacheDir == null)
                ? (Files.isDirectory(normalizedInput)
                    ? normalizedInput.resolve("cache")
                    : normalizedInput.getParent().resolve("cache"))
                : cacheDir.toAbsolutePath().normalize();
        Files.createDirectories(resolvedCacheDir);

        if (Files.isDirectory(normalizedInput)) {
            List<Path> simFiles = listSimulationFiles(normalizedInput);
            for (Path simFile : simFiles) {
                buildCache(simFile, resolvedCacheDir, dS);
            }
        } else {
            buildCache(normalizedInput, resolvedCacheDir, dS);
        }
    }

    private static List<Path> listSimulationFiles(Path directory) throws IOException {
        List<Path> files = new ArrayList<>();
        try (var paths = Files.list(directory)) {
            paths.filter(Files::isRegularFile)
                    .filter(path -> path.getFileName().toString().startsWith("sim_"))
                    .filter(path -> path.getFileName().toString().endsWith(".txt"))
                    .sorted()
                    .forEach(files::add);
        }
        return files;
    }

    private static void buildCache(Path simFile, Path cacheDir, double dS) throws IOException {
        AnalysisResult result = analyzeSimulationFile(simFile, dS);
        Path cachePath = cachePathFor(simFile, cacheDir);
        writeCache(result, cachePath, dS);
        System.out.println("Cached analysis: " + simFile.getFileName() + " -> " + cachePath);
    }

    public static Path cachePathFor(Path simFile, Path cacheDir) {
        String filename = simFile.getFileName().toString();
        String stem = filename.endsWith(".txt") ? filename.substring(0, filename.length() - 4) : filename;
        return cacheDir.resolve(stem + CACHE_SUFFIX);
    }

    private static AnalysisResult analyzeSimulationFile(Path simFile, double dS) throws IOException {
        Map<String, Double> metadata = new LinkedHashMap<>();
        List<Double> eventTimes = new ArrayList<>();
        List<Double> fuTimes = new ArrayList<>();
        List<Double> fuValues = new ArrayList<>();
        List<Double> fallbackCfcTimes = new ArrayList<>();
        List<Double> fallbackCfcValues = new ArrayList<>();

        long sourceSize = Files.size(simFile);
        FileTime sourceMtime = Files.getLastModifiedTime(simFile);

        SnapshotFrame previousFrame = null;
        ShellGeometry shellGeometry = null;
        RunningProfiles runningProfiles = null;
        int fallbackCumulative = 0;

        try (BufferedReader reader = Files.newBufferedReader(simFile, StandardCharsets.UTF_8)) {
            String line;
            int N = 0;
            while ((line = reader.readLine()) != null) {
                String trimmed = line.trim();
                if (trimmed.isEmpty()) {
                    continue;
                }

                if (trimmed.startsWith("#")) {
                    parseMetadataLine(trimmed, metadata);
                    continue;
                }

                if (trimmed.startsWith("S ")) {
                    if (shellGeometry == null) {
                        N = metadata.get("N").intValue();
                        shellGeometry = buildShellGeometry(metadata, dS);
                        runningProfiles = new RunningProfiles(shellGeometry.centers.length);
                    }

                    SnapshotFrame currentFrame = parseSnapshot(trimmed, reader, N);
                    fuTimes.add(currentFrame.time);
                    fuValues.add(currentFrame.usedCount / (double) N);

                    if (previousFrame == null) {
                        fallbackCfcTimes.add(currentFrame.time);
                        fallbackCfcValues.add(0.0);
                    } else {
                        double dt = currentFrame.time - previousFrame.time;
                        if (dt > 0.0) {
                            accumulateRadialProfiles(previousFrame, dt, shellGeometry, runningProfiles);
                        }
                        fallbackCumulative += countFreshToUsedTransitions(previousFrame.used, currentFrame.used);
                        fallbackCfcTimes.add(currentFrame.time);
                        fallbackCfcValues.add((double) fallbackCumulative);
                    }

                    previousFrame = currentFrame;
                    continue;
                }

                if (trimmed.startsWith("E ")) {
                    String[] parts = trimmed.split("\\s+");
                    eventTimes.add(Double.parseDouble(parts[1]));
                }
            }
        }

        if (shellGeometry == null || runningProfiles == null) {
            throw new IOException("No snapshots found in " + simFile);
        }

        double tFinal = metadata.getOrDefault("t_final", fuTimes.isEmpty() ? 0.0 : fuTimes.get(fuTimes.size() - 1));
        double[] cfcTimes;
        double[] cfcValues;
        if (!eventTimes.isEmpty()) {
            cfcTimes = buildExactCfcTimes(eventTimes, tFinal);
            cfcValues = buildExactCfcValues(eventTimes, tFinal);
        } else {
            cfcTimes = toDoubleArray(fallbackCfcTimes);
            cfcValues = toDoubleArray(fallbackCfcValues);
        }

        LinearFit fit = linearFit(cfcTimes, cfcValues);

        double[] rhoValues = new double[shellGeometry.centers.length];
        double[] vValues = new double[shellGeometry.centers.length];
        double[] jValues = new double[shellGeometry.centers.length];

        if (runningProfiles.totalTime > 0.0) {
            for (int i = 0; i < shellGeometry.centers.length; i++) {
                rhoValues[i] = runningProfiles.countTime[i] / (runningProfiles.totalTime * shellGeometry.areas[i]);
                if (runningProfiles.particleTime[i] > 1e-12) {
                    vValues[i] = runningProfiles.velocityTimeSum[i] / runningProfiles.particleTime[i];
                } else {
                    vValues[i] = 0.0;
                }
                jValues[i] = rhoValues[i] * Math.abs(vValues[i]);
            }
        }

        return new AnalysisResult(
                simFile.toAbsolutePath().normalize(),
                sourceSize,
                sourceMtime,
                metadata,
                cfcTimes,
                cfcValues,
                fit,
                toDoubleArray(fuTimes),
                toDoubleArray(fuValues),
                shellGeometry.centers,
                rhoValues,
                vValues,
                jValues
        );
    }

    private static void parseMetadataLine(String line, Map<String, Double> metadata) {
        String cleaned = line.replaceFirst("^#\\s*", "");
        for (String part : cleaned.split("\\s+")) {
            int eq = part.indexOf('=');
            if (eq <= 0 || eq == part.length() - 1) {
                continue;
            }
            String key = part.substring(0, eq);
            String value = part.substring(eq + 1);
            try {
                metadata.put(key, Double.parseDouble(value));
            } catch (NumberFormatException ignored) {
                // Ignore non-numeric metadata tokens.
            }
        }
    }

    private static SnapshotFrame parseSnapshot(String snapshotHeader, BufferedReader reader, int N) throws IOException {
        String[] headerParts = snapshotHeader.split("\\s+");
        double time = Double.parseDouble(headerParts[1]);

        double[] x = new double[N];
        double[] y = new double[N];
        double[] vx = new double[N];
        double[] vy = new double[N];
        boolean[] used = new boolean[N];
        int usedCount = 0;

        for (int i = 0; i < N; i++) {
            String line = reader.readLine();
            if (line == null) {
                throw new IOException("Unexpected EOF while reading snapshot at t=" + time);
            }
            String[] parts = line.trim().split("\\s+");
            int pid = Integer.parseInt(parts[0]);
            x[pid] = Double.parseDouble(parts[1]);
            y[pid] = Double.parseDouble(parts[2]);
            vx[pid] = Double.parseDouble(parts[3]);
            vy[pid] = Double.parseDouble(parts[4]);
            used[pid] = "U".equals(parts[5]);
            if (used[pid]) {
                usedCount++;
            }
        }

        return new SnapshotFrame(time, x, y, vx, vy, used, usedCount);
    }

    private static int countFreshToUsedTransitions(boolean[] previousUsed, boolean[] currentUsed) {
        int transitions = 0;
        for (int i = 0; i < previousUsed.length; i++) {
            if (!previousUsed[i] && currentUsed[i]) {
                transitions++;
            }
        }
        return transitions;
    }

    private static ShellGeometry buildShellGeometry(Map<String, Double> metadata, double dS) {
        double enclosureRadius = metadata.get("R_enclosure");
        double obstacleRadius = metadata.get("r0");
        double particleRadius = metadata.get("r");
        double sMin = obstacleRadius + particleRadius;
        double sMax = enclosureRadius - particleRadius;

        List<Double> edgeList = new ArrayList<>();
        for (double value = sMin; value <= sMax + 1e-9; value += dS) {
            edgeList.add(value);
        }
        if (edgeList.get(edgeList.size() - 1) < sMax - 1e-9) {
            edgeList.add(sMax);
        } else {
            edgeList.set(edgeList.size() - 1, sMax);
        }

        double[] edges = toDoubleArray(edgeList);
        double[] centers = new double[edges.length - 1];
        double[] areas = new double[edges.length - 1];
        for (int i = 0; i < centers.length; i++) {
            centers[i] = 0.5 * (edges[i] + edges[i + 1]);
            areas[i] = Math.PI * (edges[i + 1] * edges[i + 1] - edges[i] * edges[i]);
        }

        return new ShellGeometry(edges, centers, areas);
    }

    private static void accumulateRadialProfiles(SnapshotFrame frame, double dt, ShellGeometry geometry,
                                                 RunningProfiles runningProfiles) {
        if (dt <= 0.0) {
            return;
        }

        double[] counts = new double[geometry.centers.length];
        double[] velocitySum = new double[geometry.centers.length];

        for (int i = 0; i < frame.x.length; i++) {
            if (frame.used[i]) {
                continue;
            }

            double x = frame.x[i];
            double y = frame.y[i];
            double vx = frame.vx[i];
            double vy = frame.vy[i];
            double radialDistance = Math.sqrt(x * x + y * y);
            double rdotv = x * vx + y * vy;

            if (rdotv >= 0.0) {
                continue;
            }
            if (radialDistance < geometry.edges[0] - 1e-12 ||
                    radialDistance > geometry.edges[geometry.edges.length - 1] + 1e-12) {
                continue;
            }

            int shellIndex = upperBound(geometry.edges, radialDistance) - 1;
            shellIndex = Math.max(0, Math.min(shellIndex, geometry.centers.length - 1));
            double vfIn = radialDistance > 1e-10 ? rdotv / radialDistance : 0.0;

            counts[shellIndex] += 1.0;
            velocitySum[shellIndex] += vfIn;
        }

        for (int shellIndex = 0; shellIndex < geometry.centers.length; shellIndex++) {
            runningProfiles.countTime[shellIndex] += dt * counts[shellIndex];
            runningProfiles.velocityTimeSum[shellIndex] += dt * velocitySum[shellIndex];
            runningProfiles.particleTime[shellIndex] += dt * counts[shellIndex];
        }
        runningProfiles.totalTime += dt;
    }

    private static int upperBound(double[] values, double target) {
        int low = 0;
        int high = values.length;
        while (low < high) {
            int mid = (low + high) >>> 1;
            if (target < values[mid]) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        return low;
    }

    private static double[] buildExactCfcTimes(List<Double> eventTimes, double tFinal) {
        double[] times = new double[eventTimes.size() + 2];
        times[0] = 0.0;
        for (int i = 0; i < eventTimes.size(); i++) {
            times[i + 1] = eventTimes.get(i);
        }
        if (eventTimes.isEmpty() || eventTimes.get(eventTimes.size() - 1) < tFinal) {
            times[times.length - 1] = tFinal;
            return times;
        }
        return Arrays.copyOf(times, times.length - 1);
    }

    private static double[] buildExactCfcValues(List<Double> eventTimes, double tFinal) {
        double[] values = new double[eventTimes.size() + 2];
        values[0] = 0.0;
        for (int i = 0; i < eventTimes.size(); i++) {
            values[i + 1] = i + 1.0;
        }
        if (eventTimes.isEmpty() || eventTimes.get(eventTimes.size() - 1) < tFinal) {
            values[values.length - 1] = eventTimes.size();
            return values;
        }
        return Arrays.copyOf(values, values.length - 1);
    }

    private static LinearFit linearFit(double[] x, double[] y) {
        if (x.length < 2 || y.length < 2 || x.length != y.length) {
            return new LinearFit(0.0, 0.0, 0.0);
        }

        double meanX = 0.0;
        double meanY = 0.0;
        for (int i = 0; i < x.length; i++) {
            meanX += x[i];
            meanY += y[i];
        }
        meanX /= x.length;
        meanY /= y.length;

        double cov = 0.0;
        double varX = 0.0;
        for (int i = 0; i < x.length; i++) {
            double dx = x[i] - meanX;
            cov += dx * (y[i] - meanY);
            varX += dx * dx;
        }

        if (varX <= 1e-18) {
            return new LinearFit(0.0, meanY, 0.0);
        }

        double slope = cov / varX;
        double intercept = meanY - slope * meanX;

        double ssRes = 0.0;
        double ssTot = 0.0;
        for (int i = 0; i < x.length; i++) {
            double predicted = slope * x[i] + intercept;
            ssRes += (y[i] - predicted) * (y[i] - predicted);
            ssTot += (y[i] - meanY) * (y[i] - meanY);
        }

        double rSquared = ssTot > 0.0 ? 1.0 - ssRes / ssTot : 0.0;
        return new LinearFit(slope, intercept, rSquared);
    }

    private static double[] toDoubleArray(List<Double> values) {
        double[] array = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            array[i] = values.get(i);
        }
        return array;
    }

    private static void writeCache(AnalysisResult result, Path cachePath, double dS) throws IOException {
        Files.createDirectories(cachePath.getParent());
        try (var outputStream = Files.newOutputStream(cachePath);
             var gzipStream = new GZIPOutputStream(outputStream);
             var writer = new PrintWriter(new BufferedWriter(new java.io.OutputStreamWriter(
                     gzipStream, StandardCharsets.UTF_8)))) {

            writer.println("{");
            writer.println("  \"cache_version\": 1,");
            writer.printf("  \"source_path\": \"%s\",%n", escapeJson(result.sourcePath.toString()));
            writer.printf("  \"source_size\": %d,%n", result.sourceSize);
            writer.printf("  \"source_mtime_ms\": %d,%n", result.sourceMtime.toMillis());
            writeNumericObject(writer, "metadata", result.metadata, "  ");
            writer.println(",");
            writer.println("  \"cfc\": {");
            writeArray(writer, "times", result.cfcTimes, "    ");
            writer.println(",");
            writeArray(writer, "values", result.cfcValues, "    ");
            writer.println(",");
            writer.printf(Locale.US, "    \"J\": %.17g,%n", result.cfcFit.slope);
            writer.printf(Locale.US, "    \"intercept\": %.17g,%n", result.cfcFit.intercept);
            writer.printf(Locale.US, "    \"r_squared\": %.17g%n", result.cfcFit.rSquared);
            writer.println("  },");
            writer.println("  \"fu\": {");
            writeArray(writer, "times", result.fuTimes, "    ");
            writer.println(",");
            writeArray(writer, "values", result.fuValues, "    ");
            writer.println();
            writer.println("  },");
            writer.println("  \"radial_profiles\": {");
            writer.printf(Locale.US, "    \"dS\": %.17g,%n", dS);
            writeArray(writer, "S_centers", result.shellCenters, "    ");
            writer.println(",");
            writeArray(writer, "rho", result.rhoValues, "    ");
            writer.println(",");
            writeArray(writer, "v", result.vValues, "    ");
            writer.println(",");
            writeArray(writer, "J", result.jValues, "    ");
            writer.println();
            writer.println("  }");
            writer.println("}");
        }
    }

    private static void writeNumericObject(PrintWriter writer, String key, Map<String, Double> values, String indent) {
        writer.printf("%s\"%s\": {%n", indent, key);
        int index = 0;
        int size = values.size();
        for (Map.Entry<String, Double> entry : values.entrySet()) {
            String suffix = index + 1 < size ? "," : "";
            writer.printf(Locale.US, "%s  \"%s\": %.17g%s%n",
                    indent, escapeJson(entry.getKey()), entry.getValue(), suffix);
            index++;
        }
        writer.printf("%s}", indent);
    }

    private static void writeArray(PrintWriter writer, String key, double[] values, String indent) {
        writer.printf("%s\"%s\": [", indent, key);
        for (int i = 0; i < values.length; i++) {
            if (i > 0) {
                writer.print(", ");
            }
            writer.printf(Locale.US, "%.17g", values[i]);
        }
        writer.print("]");
    }

    private static String escapeJson(String value) {
        return value.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
