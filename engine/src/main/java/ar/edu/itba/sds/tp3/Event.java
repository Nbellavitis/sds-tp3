package ar.edu.itba.sds.tp3;


import ar.edu.itba.sds.tp3.*;
/**
 * Represents a future collision event in the Event-Driven Molecular Dynamics simulation.
 * 
 * Event types:
 * - PARTICLE_PARTICLE: collision between two particles (a and b both non-null)
 * - PARTICLE_OUTER_WALL: particle a hits the outer circular enclosure (b is null)
 * - PARTICLE_OBSTACLE: particle a hits the central fixed obstacle (b is null)
 * - SNAPSHOT: periodic state dump for output (a and b both null)
 *
 * Invalidation: each event stores the collision counts of the involved particles
 * at creation time. If the current count differs, the event is stale.
 */
public class Event implements Comparable<Event> {

    public enum Type {
        PARTICLE_PARTICLE,
        PARTICLE_OUTER_WALL,
        PARTICLE_OBSTACLE,
        SNAPSHOT
    }

    private final double time;          // absolute time of this event [s]
    private final Type type;
    private final Particle a;           // first particle (null for SNAPSHOT)
    private final Particle b;           // second particle (only for PARTICLE_PARTICLE)
    private final int countA;           // collision count of a at creation
    private final int countB;           // collision count of b at creation

    // ── Particle–Particle ────────────────────────────────────────────────
    public static Event particleParticle(double time, Particle a, Particle b) {
        return new Event(time, Type.PARTICLE_PARTICLE, a, b,
                a.getCollisionCount(), b.getCollisionCount());
    }

    // ── Particle–Outer Wall ──────────────────────────────────────────────
    public static Event particleOuterWall(double time, Particle a) {
        return new Event(time, Type.PARTICLE_OUTER_WALL, a, null,
                a.getCollisionCount(), -1);
    }

    // ── Particle–Obstacle ────────────────────────────────────────────────
    public static Event particleObstacle(double time, Particle a) {
        return new Event(time, Type.PARTICLE_OBSTACLE, a, null,
                a.getCollisionCount(), -1);
    }

    // ── Snapshot (periodic output) ───────────────────────────────────────
    public static Event snapshot(double time) {
        return new Event(time, Type.SNAPSHOT, null, null, -1, -1);
    }

    private Event(double time, Type type, Particle a, Particle b,
                  int countA, int countB) {
        this.time = time;
        this.type = type;
        this.a = a;
        this.b = b;
        this.countA = countA;
        this.countB = countB;
    }

    /**
     * An event is valid if no involved particle has had a collision
     * since this event was created.
     */
    public boolean isValid() {
        if (type == Type.SNAPSHOT) return true;
        if (a != null && a.getCollisionCount() != countA) return false;
        if (b != null && b.getCollisionCount() != countB) return false;
        return true;
    }

    @Override
    public int compareTo(Event other) {
        return Double.compare(this.time, other.time);
    }

    // ── Getters ──────────────────────────────────────────────────────────
    public double getTime() { return time; }
    public Type getType() { return type; }
    public Particle getA() { return a; }
    public Particle getB() { return b; }
}
