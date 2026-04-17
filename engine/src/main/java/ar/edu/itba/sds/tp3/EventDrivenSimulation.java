package ar.edu.itba.sds.tp3;

import java.io.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * Event-Driven Molecular Dynamics Simulation
 * System 1: Scanning rate in a circular enclosure with a fixed central obstacle.
 *
 * Enclosure: circle of diameter L=80m (radius=40m), center at (0,0).
 * Obstacle:  fixed disc at origin, radius r0=1m, infinite mass.
 * Particles: N particles, radius r=1m, mass m=1kg, initial speed v0=1 m/s.
 *
 * Output: snapshots at t=0, every configurable number of collisions, and at the end.
 */
public class EventDrivenSimulation {

    // ── Physical parameters ──────────────────────────────────────────────
    private static final double L = 80.0;                    // enclosure diameter [m]
    private static final double ENCLOSURE_RADIUS = L / 2.0;  // 40 m
    private static final double OBSTACLE_RADIUS = 1.0;       // r0 = 1 m
    private static final double PARTICLE_RADIUS = 1.0;       // r = 1 m
    private static final double PARTICLE_MASS = 1.0;         // m = 1 kg
    private static final double V0 = 1.0;                    // initial speed [m/s]
    private static final double DEFAULT_T_FINAL = 5.0;       // default simulation end time [s]
    private static final double DEFAULT_TIMING_WINDOW_1_1 = 5.0; // default simulated seconds measured for inciso 1.1
    private static final int DEFAULT_SNAPSHOT_EVERY_EVENTS = 10; // fallback keep one snapshot every N processed collisions
    private static final long TARGET_MAX_SNAPSHOT_FILE_BYTES = 100L * 1024L * 1024L;
    private static final Map<Integer, Integer> SNAPSHOT_EVERY_EVENTS_BY_N = buildSnapshotEveryEventsByN();

    // ── Data ─────────────────────────────────────────────────────────────
    private final int N;
    private final Particle[] particles;
    private final PriorityQueue<Event> pq;
    private double currentTime;
    private final Random random;
    private final double tFinal;
    private final double timingWindow;
    private final boolean writeOutput;
    private final String snapshotOutputFilePath;
    private final String transitionOutputFilePath;
    private final int snapshotEveryEvents;

    // ── Metrics ──────────────────────────────────────────────────────────
    private int totalCollisions;
    private double lastSnapshotTime;

    public static final class RunTimings {
        private final long totalElapsedNs;
        private final long elapsedToTimingWindowNs;
        private final double timingWindow;

        private RunTimings(long totalElapsedNs, long elapsedToTimingWindowNs, double timingWindow) {
            this.totalElapsedNs = totalElapsedNs;
            this.elapsedToTimingWindowNs = elapsedToTimingWindowNs;
            this.timingWindow = timingWindow;
        }

        public long getTotalElapsedNs() {
            return totalElapsedNs;
        }

        public long getElapsedToTimingWindowNs() {
            return elapsedToTimingWindowNs;
        }

        public double getTimingWindow() {
            return timingWindow;
        }
    }

    public EventDrivenSimulation(int N, long seed, double tFinal, double timingWindow,
                                 boolean writeOutput, int snapshotEveryEvents) {
        this.N = N;
        this.particles = new Particle[N];
        this.pq = new PriorityQueue<>();
        this.currentTime = 0.0;
        this.random = new Random(seed);
        this.tFinal = tFinal;
        this.timingWindow = timingWindow;
        this.writeOutput = writeOutput;
        this.snapshotEveryEvents = snapshotEveryEvents;
        this.totalCollisions = 0;
        this.lastSnapshotTime = Double.NaN;

        if (snapshotEveryEvents <= 0) {
            throw new IllegalArgumentException("snapshot_every_events must be >= 1");
        }

        if (writeOutput) {
            String timestamp = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss"));
            String suffix = N + "N_" + timestamp + "_s" + seed + ".txt";
            this.snapshotOutputFilePath = "data/sim_" + suffix;
            this.transitionOutputFilePath = "data/events_" + suffix;
        } else {
            this.snapshotOutputFilePath = null;
            this.transitionOutputFilePath = null;
        }

        initializeParticles();
    }

    private static Map<Integer, Integer> buildSnapshotEveryEventsByN() {
        Map<Integer, Integer> values = new LinkedHashMap<>();
        // Manual per-N overrides. Edit freely if a specific N needs denser or sparser snapshots.
        values.put(100, 10);
        values.put(150, 15);
        values.put(200, 20);
        values.put(250, 25);
        values.put(300, 30);
        values.put(350, 35);
        values.put(400, 40);
        values.put(450, 45);
        values.put(500, 50);
        values.put(550, 60);
        values.put(600, 70);
        values.put(650, 80);
        values.put(700, 90);
        values.put(750, 100);
        values.put(800, 120);
        return Collections.unmodifiableMap(values);
    }

    private static int resolveDefaultSnapshotEveryEvents(int n) {
        Integer manual = SNAPSHOT_EVERY_EVENTS_BY_N.get(n);
        if (manual != null && manual > 0) {
            return manual;
        }
        return DEFAULT_SNAPSHOT_EVERY_EVENTS;
    }

    // ── Particle initialization ──────────────────────────────────────────
    private void initializeParticles() {
        int placed = 0;
        int id = 0;

        while (placed < N) {
            // Random position inside the circular enclosure, outside the obstacle
            double angle = random.nextDouble() * 2.0 * Math.PI;
            double minR = OBSTACLE_RADIUS + PARTICLE_RADIUS;
            double maxR = ENCLOSURE_RADIUS - PARTICLE_RADIUS;
            // Uniform distribution in area: r = sqrt(U * (maxR^2 - minR^2) + minR^2)
            double r = Math.sqrt(random.nextDouble() * (maxR * maxR - minR * minR) + minR * minR);

            double x = r * Math.cos(angle);
            double y = r * Math.sin(angle);

            // Check overlap with already placed particles
            boolean overlap = false;
            for (int i = 0; i < placed; i++) {
                double dx = x - particles[i].getX();
                double dy = y - particles[i].getY();
                double dist = Math.sqrt(dx * dx + dy * dy);
                if (dist < 2.0 * PARTICLE_RADIUS) {
                    overlap = true;
                    break;
                }
            }
            if (overlap) continue;

            // Random velocity direction, uniform in [0, 2pi)
            double velAngle = random.nextDouble() * 2.0 * Math.PI;
            double vx = V0 * Math.cos(velAngle);
            double vy = V0 * Math.sin(velAngle);

            particles[placed] = new Particle(id, x, y, vx, vy, PARTICLE_RADIUS, PARTICLE_MASS);
            placed++;
            id++;
        }
    }

    // ── Predict all events for particle i ────────────────────────────────
    private void predictEvents(int i) {
        Particle pi = particles[i];

        // Particle-Particle collisions
        for (int j = 0; j < N; j++) {
            if (j == i) continue;
            double dt = pi.timeToCollide(particles[j], currentTime);
            if (Double.isFinite(dt) && currentTime + dt < tFinal) {
                pq.add(Event.particleParticle(currentTime + dt, pi, particles[j]));
            }
        }

        // Outer wall collision
        double dtWall = pi.timeToOuterWall(ENCLOSURE_RADIUS, currentTime);
        if (Double.isFinite(dtWall) && currentTime + dtWall < tFinal) {
            pq.add(Event.particleOuterWall(currentTime + dtWall, pi));
        }

        // Obstacle collision
        double dtObs = pi.timeToObstacle(OBSTACLE_RADIUS, currentTime);
        if (Double.isFinite(dtObs) && currentTime + dtObs < tFinal) {
            pq.add(Event.particleObstacle(currentTime + dtObs, pi));
        }
    }

    // ── Run the simulation ───────────────────────────────────────────────
    public RunTimings run() {
        if (writeOutput) {
            new File("data").mkdirs();
        }

        long startWall = System.nanoTime();
        double effectiveTimingWindow = Math.min(timingWindow, tFinal);
        long elapsedToTimingWindowNs = -1L;

        // Initialize PQ with all events
        for (int i = 0; i < N; i++) {
            predictEvents(i);
        }

        PrintWriter snapshotWriter = null;
        PrintWriter transitionWriter = null;
        try {
            if (writeOutput) {
                snapshotWriter = new PrintWriter(new BufferedWriter(new FileWriter(snapshotOutputFilePath)));
                transitionWriter = new PrintWriter(new BufferedWriter(new FileWriter(transitionOutputFilePath)));
            }

            if (snapshotWriter != null && transitionWriter != null) {
                snapshotWriter.println("# N=" + N + " L=" + L + " R_enclosure=" + ENCLOSURE_RADIUS
                        + " r0=" + OBSTACLE_RADIUS + " r=" + PARTICLE_RADIUS
                        + " m=" + PARTICLE_MASS + " v0=" + V0 + " t_final=" + tFinal
                        + " snapshot_every_events=" + snapshotEveryEvents);
                snapshotWriter.println("# transition_log_file=" + transitionOutputFilePath);
                snapshotWriter.println("# FORMAT: SNAPSHOT lines start with 'S', followed by time,");
                snapshotWriter.println("# then N lines of: id x y vx vy state(F/U)");

                transitionWriter.println("# N=" + N + " L=" + L + " R_enclosure=" + ENCLOSURE_RADIUS
                        + " r0=" + OBSTACLE_RADIUS + " r=" + PARTICLE_RADIUS
                        + " m=" + PARTICLE_MASS + " v0=" + V0 + " t_final=" + tFinal
                        + " snapshot_every_events=" + snapshotEveryEvents);
                transitionWriter.println("# snapshot_file=" + snapshotOutputFilePath);
                transitionWriter.println("# FORMAT: T time id from_state to_state");

                writeSnapshot(snapshotWriter);
            }

            // Main event loop
            while (!pq.isEmpty()) {
                Event event = pq.poll();

                // Skip invalid (stale) events
                if (!event.isValid()) continue;

                if (elapsedToTimingWindowNs < 0
                        && currentTime < effectiveTimingWindow
                        && event.getTime() > effectiveTimingWindow) {
                    currentTime = effectiveTimingWindow;
                    elapsedToTimingWindowNs = System.nanoTime() - startWall;
                }

                // Event time beyond simulation end
                if (event.getTime() > tFinal) break;

                currentTime = event.getTime();
                if (elapsedToTimingWindowNs < 0 && currentTime >= effectiveTimingWindow) {
                    elapsedToTimingWindowNs = System.nanoTime() - startWall;
                }

                // Process event
                boolean processedCollision = false;
                switch (event.getType()) {
                    case PARTICLE_PARTICLE -> {
                        Particle a = event.getA();
                        Particle b = event.getB();
                        a.resolveCollision(b, currentTime);
                        totalCollisions++;
                        processedCollision = true;

                        // Re-predict events for both particles
                        int ia = findIndex(a);
                        int ib = findIndex(b);
                        predictEvents(ia);
                        predictEvents(ib);
                    }
                    case PARTICLE_OUTER_WALL -> {
                        Particle a = event.getA();
                        boolean wasUsedToFresh = a.resolveOuterWallCollision(currentTime);
                        totalCollisions++;
                        processedCollision = true;

                        if (wasUsedToFresh && transitionWriter != null) {
                            transitionWriter.printf("T %.6e %d U F%n", currentTime, a.getId());
                        }

                        int ia = findIndex(a);
                        predictEvents(ia);
                    }
                    case PARTICLE_OBSTACLE -> {
                        Particle a = event.getA();
                        boolean wasFreshToUsed = a.resolveObstacleCollision(currentTime);
                        totalCollisions++;
                        processedCollision = true;

                        if (wasFreshToUsed && transitionWriter != null) {
                            transitionWriter.printf("T %.6e %d F U%n", currentTime, a.getId());
                        }

                        int ia = findIndex(a);
                        predictEvents(ia);
                    }
                    case SNAPSHOT -> {
                        // Snapshot events are no longer scheduled. Keep the case
                        // to avoid surprising future callers if they are reintroduced.
                    }
                }

                if (snapshotWriter != null && processedCollision && totalCollisions % snapshotEveryEvents == 0) {
                    writeSnapshot(snapshotWriter);
                }
            }

            // Final snapshot
            if (currentTime < tFinal) {
                currentTime = tFinal;
            }

            if (elapsedToTimingWindowNs < 0 && currentTime >= effectiveTimingWindow) {
                elapsedToTimingWindowNs = System.nanoTime() - startWall;
            }

            if (snapshotWriter != null && (Double.isNaN(lastSnapshotTime) || Math.abs(lastSnapshotTime - currentTime) > 1e-12)) {
                writeSnapshot(snapshotWriter);
            }

        } catch (IOException e) {
            System.err.println("Error writing output: " + e.getMessage());
            e.printStackTrace();
        } finally {
            if (snapshotWriter != null) {
                snapshotWriter.close();
            }
            if (transitionWriter != null) {
                transitionWriter.close();
            }
        }

        long endWall = System.nanoTime();
        long totalElapsedNs = endWall - startWall;
        if (elapsedToTimingWindowNs < 0) {
            elapsedToTimingWindowNs = totalElapsedNs;
        }
        return new RunTimings(totalElapsedNs, elapsedToTimingWindowNs, effectiveTimingWindow);
    }

    // ── Write snapshot of all particles ──────────────────────────────────
    private void writeSnapshot(PrintWriter writer) {
        writer.printf("S %.6e%n", currentTime);
        for (Particle p : particles) {
            writer.printf("%d %.6e %.6e %.6e %.6e %s%n",
                    p.getId(), p.getX(currentTime), p.getY(currentTime), p.getVx(), p.getVy(),
                    p.getState() == Particle.State.FRESH ? "F" : "U");
        }
        lastSnapshotTime = currentTime;
    }

    // ── Find index of particle in array ──────────────────────────────────
    private int findIndex(Particle p) {
        // Since particle IDs are assigned 0..N-1 and array is indexed the same
        return p.getId();
    }

    // ── Main entry point ─────────────────────────────────────────────────
    public static void main(String[] args) {
        int N = 100;          // default
        long seed = -1;       // -1 means random seed
        int runs = 1;         // number of independent realizations
        double tFinal = DEFAULT_T_FINAL;
        double timingWindow = DEFAULT_TIMING_WINDOW_1_1;
        boolean writeOutput = true;
        Integer snapshotEveryEventsOverride = null;

        // Parse command-line arguments
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-N" -> N = Integer.parseInt(args[++i]);
                case "-seed" -> seed = Long.parseLong(args[++i]);
                case "-runs" -> runs = Integer.parseInt(args[++i]);
                case "-t", "-tf", "-t_final", "-tfinal" -> tFinal = Double.parseDouble(args[++i]);
                case "--timing-window" -> timingWindow = Double.parseDouble(args[++i]);
                case "-snapshot_every_events", "-snapshot-every-events", "--snapshot-every-events" ->
                        snapshotEveryEventsOverride = Integer.parseInt(args[++i]);
                case "--no-output", "--timing-only" -> writeOutput = false;
            }
        }

        int snapshotEveryEvents = writeOutput
                ? resolveDefaultSnapshotEveryEvents(N)
                : DEFAULT_SNAPSHOT_EVERY_EVENTS;
        if (snapshotEveryEventsOverride != null) {
            snapshotEveryEvents = snapshotEveryEventsOverride;
        }

        if (!Double.isFinite(tFinal) || tFinal <= 0.0) {
            throw new IllegalArgumentException("t_final must be a positive finite value");
        }
        if (!Double.isFinite(timingWindow) || timingWindow <= 0.0) {
            throw new IllegalArgumentException("timing_window must be a positive finite value");
        }
        if (snapshotEveryEvents <= 0) {
            throw new IllegalArgumentException("snapshot_every_events must be >= 1");
        }

        System.out.println("=== Event-Driven Molecular Dynamics ===");
        System.out.println("System 1: Scanning rate in circular enclosure");
        System.out.println("N = " + N + ", L = " + L + " m, t_final = " + tFinal + " s");
        System.out.println("Enclosure radius = " + ENCLOSURE_RADIUS + " m");
        System.out.println("Obstacle radius  = " + OBSTACLE_RADIUS + " m");
        System.out.println("Particle radius  = " + PARTICLE_RADIUS + " m, mass = " + PARTICLE_MASS + " kg");
        System.out.println("Initial speed    = " + V0 + " m/s");
        System.out.println("Runs = " + runs);
        if (writeOutput) {
            System.out.println("Snapshot every   = " + snapshotEveryEvents + " collisions");
            if (snapshotEveryEventsOverride != null) {
                System.out.println("Snapshot policy  = CLI override");
            } else if (SNAPSHOT_EVERY_EVENTS_BY_N.containsKey(N)) {
                System.out.println("Snapshot policy  = manual per-N map");
            } else {
                System.out.println("Snapshot policy  = default fallback");
            }
        }
        System.out.println();

        for (int run = 0; run < runs; run++) {
            long actualSeed = (seed < 0) ? System.nanoTime() : seed + run;
            System.out.println("--- Run " + (run + 1) + "/" + runs + " (seed=" + actualSeed + ") ---");

            EventDrivenSimulation sim = new EventDrivenSimulation(
                    N,
                    actualSeed,
                    tFinal,
                    timingWindow,
                    writeOutput,
                    snapshotEveryEvents
            );
            RunTimings timings = sim.run();
            double totalMs = timings.getTotalElapsedNs() / 1e6;
            double timingWindowMs = timings.getElapsedToTimingWindowNs() / 1e6;

            System.out.printf("  Total runtime: %.2f ms (%.4f s)%n", totalMs, totalMs / 1000.0);
            System.out.printf("  Inciso 1.1 runtime (first %.2f simulated seconds): %.2f ms (%.4f s)%n",
                    timings.getTimingWindow(), timingWindowMs, timingWindowMs / 1000.0);
            System.out.println("  Total collisions: " + sim.totalCollisions);
            if (writeOutput) {
                System.out.println("  Output: " + sim.snapshotOutputFilePath);
                System.out.println("  Events: " + sim.transitionOutputFilePath);
                try {
                    long snapshotBytes = Files.size(Path.of(sim.snapshotOutputFilePath));
                    double snapshotMiB = snapshotBytes / (1024.0 * 1024.0);
                    System.out.printf("  Snapshot file size: %.2f MiB%n", snapshotMiB);
                    if (snapshotBytes > TARGET_MAX_SNAPSHOT_FILE_BYTES) {
                        System.out.printf(
                                "  WARNING: snapshot file exceeded %.1f MiB; raise snapshot_every_events for N=%d.%n",
                                TARGET_MAX_SNAPSHOT_FILE_BYTES / (1024.0 * 1024.0),
                                N
                        );
                    }
                } catch (IOException e) {
                    System.out.println("  WARNING: could not measure snapshot file size: " + e.getMessage());
                }
            } else {
                System.out.println("  Output: disabled (--no-output)");
            }
            System.out.println();
        }
    }
}
