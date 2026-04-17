package ar.edu.itba.sds.tp3;

/**
 * Particle for Event-Driven Molecular Dynamics in a circular enclosure.
 * Enclosure: circle of diameter L centered at origin.
 * Obstacle: fixed disc at origin with radius r0.
 * States: FRESH (green) / USED (violet).
 */
public class Particle {

    public enum State { FRESH, USED }

    private final int id;
    private double x, y;       // position [m]
    private double vx, vy;     // velocity [m/s]
    private final double radius;     // particle radius [m]
    private final double mass;       // particle mass [kg]
    private double lastUpdateTime;   // absolute time of stored (x, y)
    private int collisionCount;      // for event invalidation
    private State state;

    public Particle(int id, double x, double y, double vx, double vy,
                    double radius, double mass) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.radius = radius;
        this.mass = mass;
        this.lastUpdateTime = 0.0;
        this.collisionCount = 0;
        this.state = State.FRESH;
    }

    // ── Advance particle to time dt ──────────────────────────────────────
    public void advance(double dt) {
        if (dt <= 0.0) {
            return;
        }
        x += vx * dt;
        y += vy * dt;
        lastUpdateTime += dt;
    }

    public void advanceTo(double time) {
        if (time <= lastUpdateTime) {
            return;
        }
        advance(time - lastUpdateTime);
    }

    // ── Collision time: Particle–Particle ────────────────────────────────
    /**
     * Time until this particle collides with another particle.
     * Returns Double.POSITIVE_INFINITY if no collision.
     * Uses the quadratic from Slide 14/20 of Teorica_3.
     */
    public double timeToCollide(Particle other, double referenceTime) {
        double dx = other.getX(referenceTime) - this.getX(referenceTime);
        double dy = other.getY(referenceTime) - this.getY(referenceTime);
        double dvx = other.vx - this.vx;
        double dvy = other.vy - this.vy;

        double dvdr = dvx * dx + dvy * dy;
        if (dvdr >= 0) return Double.POSITIVE_INFINITY; // particles diverging

        double dvdv = dvx * dvx + dvy * dvy;
        double drdr = dx * dx + dy * dy;
        double sigma = this.radius + other.radius;

        double d = dvdr * dvdr - dvdv * (drdr - sigma * sigma);
        if (d < 0) return Double.POSITIVE_INFINITY; // no collision

        double tc = -(dvdr + Math.sqrt(d)) / dvdv;
        if (tc < 1e-12) return Double.POSITIVE_INFINITY; // collision in the past or now
        return tc;
    }

    // ── Collision time: Outer circular wall ──────────────────────────────
    /**
     * Time until this particle hits the outer circular wall.
     * Enclosure is a circle of radius R_outer centered at origin.
     * Particle center must not go beyond R_outer - radius.
     *
     * Solve: |r + v*t|^2 = (R_outer - radius)^2
     * (v·v)t^2 + 2(r·v)t + (r·r - R_eff^2) = 0
     * We want the POSITIVE root (the larger one if both are positive,
     * since the particle is inside).
     */
    public double timeToOuterWall(double enclosureRadius, double referenceTime) {
        double Reff = enclosureRadius - this.radius;
        double px = getX(referenceTime);
        double py = getY(referenceTime);

        double a = vx * vx + vy * vy;
        if (a < 1e-15) return Double.POSITIVE_INFINITY; // stationary

        double b = px * vx + py * vy; // half of the linear coefficient
        double c = px * px + py * py - Reff * Reff;

        double disc = b * b - a * c;
        if (disc < 0) return Double.POSITIVE_INFINITY;

        // Two roots: t = (-b ± sqrt(disc)) / a
        double sqrtDisc = Math.sqrt(disc);
        double t1 = (-b - sqrtDisc) / a;
        double t2 = (-b + sqrtDisc) / a;

        // Particle is inside the enclosure, so c <= 0.
        // We need the positive root. t2 >= t1.
        // If particle is moving outward (b > 0), both roots: t1 < 0 < t2 typically.
        // We want the first future positive crossing.
        double tc = Double.POSITIVE_INFINITY;
        if (t2 > 1e-12) tc = t2;
        if (t1 > 1e-12 && t1 < tc) tc = t1;

        return tc;
    }

    // ── Collision time: Inner obstacle (fixed disc at origin) ────────────
    /**
     * Time until this particle hits the fixed obstacle at the origin.
     * Obstacle has radius r0, particle has radius r.
     * Contact when |r + v*t| = r0 + radius.
     *
     * Solve: (v·v)t^2 + 2(r·v)t + (r·r - sigma^2) = 0
     * We want the SMALLEST POSITIVE root (particle approaches obstacle from outside).
     */
    public double timeToObstacle(double obstacleRadius, double referenceTime) {
        double sigma = obstacleRadius + this.radius;
        double px = getX(referenceTime);
        double py = getY(referenceTime);

        double a = vx * vx + vy * vy;
        if (a < 1e-15) return Double.POSITIVE_INFINITY;

        double b = px * vx + py * vy; // half of linear coefficient (= r·v)
        double c = px * px + py * py - sigma * sigma;

        // c > 0 means particle center is outside contact distance (normal case)
        // c <= 0 would be overlap or touching
        double disc = b * b - a * c;
        if (disc < 0) return Double.POSITIVE_INFINITY; // trajectory misses obstacle

        double sqrtDisc = Math.sqrt(disc);
        double t1 = (-b - sqrtDisc) / a;
        double t2 = (-b + sqrtDisc) / a;

        // We want the smallest positive root.
        // t1 <= t2 always. For a particle outside (c > 0):
        //   - If both roots positive: t1 is the entry time
        //   - If t1 < 0 < t2: particle starts inside (shouldn't happen)
        // The particle must also be approaching: we check via positive root existence.
        if (t1 > 1e-12) return t1;
        if (t2 > 1e-12) return t2;
        return Double.POSITIVE_INFINITY;
    }

    // ── Resolve elastic collision: Particle–Particle ─────────────────────
    /**
     * Elastic collision (cn=1, ct=1) using impulse J from Slide 36-37.
     * J = 2 * m1 * m2 * (dv · dr) / (sigma * (m1 + m2))
     * Updates velocities of both particles.
     */
    public void resolveCollision(Particle other, double eventTime) {
        advanceTo(eventTime);
        other.advanceTo(eventTime);

        double dx = other.x - this.x;
        double dy = other.y - this.y;
        double dvx = other.vx - this.vx;
        double dvy = other.vy - this.vy;

        double dvdr = dvx * dx + dvy * dy;
        double sigma = this.radius + other.radius;

        double J = (2.0 * this.mass * other.mass * dvdr) / (sigma * (this.mass + other.mass));

        double Jx = J * dx / sigma;
        double Jy = J * dy / sigma;

        this.vx += Jx / this.mass;
        this.vy += Jy / this.mass;
        other.vx -= Jx / other.mass;
        other.vy -= Jy / other.mass;

        this.collisionCount++;
        other.collisionCount++;
    }

    // ── Resolve collision with outer circular wall ───────────────────────
    /**
     * Specular reflection off the inner surface of the circular enclosure.
     * Normal vector points inward: n = -r/|r|
     * v' = v - 2*(v · n_out)*n_out  where n_out = r/|r|
     */
    public boolean resolveOuterWallCollision(double eventTime) {
        advanceTo(eventTime);
        double dist = Math.sqrt(x * x + y * y);
        if (dist < 1e-15) return false;

        // Outward normal at contact point
        double nx = x / dist;
        double ny = y / dist;

        double vn = vx * nx + vy * ny;
        vx -= 2.0 * vn * nx;
        vy -= 2.0 * vn * ny;

        this.collisionCount++;
        // Hitting outer wall -> becomes FRESH
        boolean wasUsed = (this.state == State.USED);
        this.state = State.FRESH;
        return wasUsed;
    }

    // ── Resolve collision with fixed obstacle at center ──────────────────
    /**
     * Specular reflection off the obstacle surface.
     * Normal vector points outward from obstacle: n = r/|r|
     * v' = v - 2*(v · n)*n
     *
     * Returns true if this was a FRESH -> USED transition (for C_fc counting).
     */
    public boolean resolveObstacleCollision(double eventTime) {
        advanceTo(eventTime);
        double dist = Math.sqrt(x * x + y * y);
        if (dist < 1e-15) return false;

        // Outward normal from obstacle (points away from origin toward particle)
        double nx = x / dist;
        double ny = y / dist;

        double vn = vx * nx + vy * ny;
        vx -= 2.0 * vn * nx;
        vy -= 2.0 * vn * ny;

        this.collisionCount++;

        // Only FRESH particles transition to USED (for C_fc counting)
        boolean wasFresh = (this.state == State.FRESH);
        this.state = State.USED; // Always USED after hitting obstacle
        return wasFresh;
    }

    // ── Getters & Setters ────────────────────────────────────────────────
    public int getId() { return id; }
    public double getX() { return x; }
    public double getY() { return y; }
    public double getX(double time) { return x + vx * (time - lastUpdateTime); }
    public double getY(double time) { return y + vy * (time - lastUpdateTime); }
    public double getVx() { return vx; }
    public double getVy() { return vy; }
    public double getRadius() { return radius; }
    public double getMass() { return mass; }
    public int getCollisionCount() { return collisionCount; }
    public State getState() { return state; }

    public void setX(double x) { this.x = x; }
    public void setY(double y) { this.y = y; }
    public void setVx(double vx) { this.vx = vx; }
    public void setVy(double vy) { this.vy = vy; }
    public void setState(State state) { this.state = state; }

    @Override
    public String toString() {
        return String.format("P[%d](%.4f,%.4f) v=(%.4f,%.4f) %s",
                id, x, y, vx, vy, state);
    }
}
