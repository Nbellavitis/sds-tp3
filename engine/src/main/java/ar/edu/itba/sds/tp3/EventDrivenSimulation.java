package ar.edu.itba.sds.tp3;

import java.io.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;

/**
 * Event-Driven Molecular Dynamics Simulation
 * System 1: Scanning rate in a circular enclosure with a fixed central obstacle.
 *
 * Enclosure: circle of diameter L=80m (radius=40m), center at (0,0).
 * Obstacle:  fixed disc at origin, radius r0=1m, infinite mass.
 * Particles: N particles, radius r=1m, mass m=1kg, initial speed v0=1 m/s.
 *
 * Output: All particle states at each collision event, written to data/ folder.
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
    private static final double TIMING_WINDOW_1_1 = 5.0;     // simulated seconds measured for inciso 1.1
    private static final int SNAPSHOT_EVERY_EVENTS = 10;     // keep one snapshot every N processed collisions

    // ── Data ─────────────────────────────────────────────────────────────
    private final int N;
    private final Particle[] particles;
    private final PriorityQueue<Event> pq;
    private double currentTime;
    private final Random random;
    private final double tFinal;
    private final String outputFilePath;

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

    public EventDrivenSimulation(int N, long seed, double tFinal) {
        this.N = N;
        this.particles = new Particle[N];
        this.pq = new PriorityQueue<>();
        this.currentTime = 0.0;
        this.random = new Random(seed);
        this.tFinal = tFinal;
        this.totalCollisions = 0;
        this.lastSnapshotTime = Double.NaN;

        // Generate unique filename with timestamp and seed
        String timestamp = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss"));
        this.outputFilePath = "data/sim_" + N + "N_" + timestamp + "_s" + seed + ".txt";

        initializeParticles();
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
            double dt = pi.timeToCollide(particles[j]);
            if (currentTime + dt < tFinal) {
                pq.add(Event.particleParticle(currentTime + dt, pi, particles[j]));
            }
        }

        // Outer wall collision
        double dtWall = pi.timeToOuterWall(ENCLOSURE_RADIUS);
        if (currentTime + dtWall < tFinal) {
            pq.add(Event.particleOuterWall(currentTime + dtWall, pi));
        }

        // Obstacle collision
        double dtObs = pi.timeToObstacle(OBSTACLE_RADIUS);
        if (currentTime + dtObs < tFinal) {
            pq.add(Event.particleObstacle(currentTime + dtObs, pi));
        }
    }

    // ── Run the simulation ───────────────────────────────────────────────
    public RunTimings run() {
        // Ensure data directory exists
        new File("data").mkdirs();

        long startWall = System.nanoTime();
        double timingWindow = Math.min(TIMING_WINDOW_1_1, tFinal);
        long elapsedToTimingWindowNs = -1L;

        // Initialize PQ with all events
        for (int i = 0; i < N; i++) {
            predictEvents(i);
        }

        try (PrintWriter writer = new PrintWriter(new BufferedWriter(
                new FileWriter(outputFilePath)))) {

            // Write metadata header
            writer.println("# N=" + N + " L=" + L + " R_enclosure=" + ENCLOSURE_RADIUS
                    + " r0=" + OBSTACLE_RADIUS + " r=" + PARTICLE_RADIUS
                    + " m=" + PARTICLE_MASS + " v0=" + V0 + " t_final=" + tFinal
                    + " snapshot_every_events=" + SNAPSHOT_EVERY_EVENTS);
            writer.println("# FORMAT: SNAPSHOT lines start with 'S', followed by time,");
            writer.println("# then N lines of: id x y vx vy state(F/U)");
            writer.println("# EVENT lines start with 'E': time id");

            // Initial state at t = 0
            writeSnapshot(writer);

            // Main event loop
            while (!pq.isEmpty()) {
                Event event = pq.poll();

                // Skip invalid (stale) events
                if (!event.isValid()) continue;

                if (elapsedToTimingWindowNs < 0 && currentTime < timingWindow && event.getTime() > timingWindow) {
                    advanceAllParticles(timingWindow - currentTime);
                    currentTime = timingWindow;
                    elapsedToTimingWindowNs = System.nanoTime() - startWall;
                }

                // Event time beyond simulation end
                if (event.getTime() > tFinal) break;

                // Advance all particles to event time
                double dt = event.getTime() - currentTime;
                advanceAllParticles(dt);
                currentTime = event.getTime();
                if (elapsedToTimingWindowNs < 0 && currentTime >= timingWindow) {
                    elapsedToTimingWindowNs = System.nanoTime() - startWall;
                }

                // Process event
                boolean processedCollision = false;
                switch (event.getType()) {
                    case PARTICLE_PARTICLE -> {
                        Particle a = event.getA();
                        Particle b = event.getB();
                        a.resolveCollision(b);
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
                        a.resolveOuterWallCollision(); // sets state to FRESH
                        totalCollisions++;
                        processedCollision = true;

                        int ia = findIndex(a);
                        predictEvents(ia);
                    }
                    case PARTICLE_OBSTACLE -> {
                        Particle a = event.getA();
                        boolean wasFreshToUsed = a.resolveObstacleCollision();
                        totalCollisions++;
                        processedCollision = true;

                        // Log F->U transition for C_fc counting
                        if (wasFreshToUsed) {
                            writer.printf("E %.6e %d%n", currentTime, a.getId());
                        }

                        int ia = findIndex(a);
                        predictEvents(ia);
                    }
                    case SNAPSHOT -> {
                        // Snapshot events are no longer scheduled. Keep the case
                        // to avoid surprising future callers if they are reintroduced.
                    }
                }

                if (processedCollision && totalCollisions % SNAPSHOT_EVERY_EVENTS == 0) {
                    writeSnapshot(writer);
                }
            }

            // Final snapshot
            if (currentTime < tFinal) {
                advanceAllParticles(tFinal - currentTime);
                currentTime = tFinal;
            }

            if (elapsedToTimingWindowNs < 0 && currentTime >= timingWindow) {
                elapsedToTimingWindowNs = System.nanoTime() - startWall;
            }

            if (Double.isNaN(lastSnapshotTime) || Math.abs(lastSnapshotTime - currentTime) > 1e-12) {
                writeSnapshot(writer);
            }

        } catch (IOException e) {
            System.err.println("Error writing output: " + e.getMessage());
            e.printStackTrace();
        }

        long endWall = System.nanoTime();
        long totalElapsedNs = endWall - startWall;
        if (elapsedToTimingWindowNs < 0) {
            elapsedToTimingWindowNs = totalElapsedNs;
        }
        return new RunTimings(totalElapsedNs, elapsedToTimingWindowNs, timingWindow);
    }

    // ── Write snapshot of all particles ──────────────────────────────────
    private void writeSnapshot(PrintWriter writer) {
        writer.printf("S %.6e%n", currentTime);
        for (Particle p : particles) {
            writer.printf("%d %.6e %.6e %.6e %.6e %s%n",
                    p.getId(), p.getX(), p.getY(), p.getVx(), p.getVy(),
                    p.getState() == Particle.State.FRESH ? "F" : "U");
        }
        lastSnapshotTime = currentTime;
    }

    private void advanceAllParticles(double dt) {
        if (dt <= 0.0) {
            return;
        }
        for (Particle p : particles) {
            p.advance(dt);
        }
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

        // Parse command-line arguments
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-N" -> N = Integer.parseInt(args[++i]);
                case "-seed" -> seed = Long.parseLong(args[++i]);
                case "-runs" -> runs = Integer.parseInt(args[++i]);
                case "-t", "-tf", "-t_final", "-tfinal" -> tFinal = Double.parseDouble(args[++i]);
            }
        }

        if (!Double.isFinite(tFinal) || tFinal <= 0.0) {
            throw new IllegalArgumentException("t_final must be a positive finite value");
        }

        System.out.println("=== Event-Driven Molecular Dynamics ===");
        System.out.println("System 1: Scanning rate in circular enclosure");
        System.out.println("N = " + N + ", L = " + L + " m, t_final = " + tFinal + " s");
        System.out.println("Enclosure radius = " + ENCLOSURE_RADIUS + " m");
        System.out.println("Obstacle radius  = " + OBSTACLE_RADIUS + " m");
        System.out.println("Particle radius  = " + PARTICLE_RADIUS + " m, mass = " + PARTICLE_MASS + " kg");
        System.out.println("Initial speed    = " + V0 + " m/s");
        System.out.println("Runs = " + runs);
        System.out.println();

        for (int run = 0; run < runs; run++) {
            long actualSeed = (seed < 0) ? System.nanoTime() : seed + run;
            System.out.println("--- Run " + (run + 1) + "/" + runs + " (seed=" + actualSeed + ") ---");

            EventDrivenSimulation sim = new EventDrivenSimulation(N, actualSeed, tFinal);
            RunTimings timings = sim.run();
            double elapsedMs = timings.getTotalElapsedNs() / 1e6;
            double timingWindowMs = timings.getElapsedToTimingWindowNs() / 1e6;

            System.out.printf("  Completed first %.2f simulated seconds in %.2f ms (%.4f s)%n",
                    timings.getTimingWindow(), timingWindowMs, timingWindowMs / 1000.0);
            System.out.printf("  Completed full run in %.2f ms (%.4f s)%n", elapsedMs, elapsedMs / 1000.0);
            System.out.println("  Total collisions: " + sim.totalCollisions);
            System.out.println("  Output: " + sim.outputFilePath);
            System.out.println();
        }
    }
}
