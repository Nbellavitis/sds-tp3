# TP3 — Dinámica Molecular Dirigida por Eventos

## Sistema 1: Scanning Rate en Recinto Cerrado con Obstáculo Fijo

Simulación de partículas en un recinto circular con un obstáculo fijo en el centro, utilizando Dinámica Molecular Dirigida por Eventos (Event-Driven Molecular Dynamics).

---

## Estructura del Proyecto

```
sds-tp3/
├── engine/                          # Motor de simulación (Java 17)
│   ├── pom.xml
│   └── src/main/java/ar/edu/itba/sds/tp3/
│       ├── Particle.java            # Partícula con colisiones y estados
│       ├── Event.java               # Eventos para la PriorityQueue
│       └── EventDrivenSimulation.java  # Motor principal (main)
│
├── graphics/                        # Post-procesamiento (Python 3)
│   ├── plot_metrics.py              # Script principal (orquestador)
│   ├── plot_execution_time.py       # Inciso 1.1: Tiempo de ejecución vs N
│   ├── plot_scanning_rate.py        # Inciso 1.2: Scanning rate <J>(N)
│   ├── plot_fraction_used.py        # Inciso 1.3: Fracción de usadas F_u(t)
│   ├── plot_radial_profiles.py      # Inciso 1.4: Perfiles radiales
│   ├── animate_system1.py           # Animación del sistema
│   └── run_batch.py                 # Runner de batch para múltiples N
│
├── data/                            # Archivos de salida de la simulación
│   └── sim_<N>N_<timestamp>.txt     # Formato: snapshots periódicos
│
└── README.md                        # Este archivo
```

---

## Requisitos

### Java
- **JDK 17** o superior
- **Maven 3.x**

### Python
- **Python 3.8+**
- **Dependencias:**
  ```bash
  pip install numpy matplotlib
  ```

---

## Compilación (Java)

```bash
cd engine
mvn compile
```

Esto compila las clases en `engine/target/classes/`.

---

## Ejecución de la Simulación

### Ejecución simple (un valor de N)

```bash
java -cp engine/target/classes ar.edu.itba.sds.tp3.EventDrivenSimulation -N 100
```

### Parámetros disponibles

| Parámetro | Descripción                     | Default |
|-----------|---------------------------------|---------|
| `-N`      | Número de partículas            | 100     |
| `-seed`   | Semilla para reproducibilidad   | random  |
| `-runs`   | Número de realizaciones         | 1       |

### Ejemplo con múltiples runs

```bash
java -cp engine/target/classes ar.edu.itba.sds.tp3.EventDrivenSimulation -N 200 -runs 5 -seed 42
```

### Ejecución batch (múltiples N)

```bash
python graphics/run_batch.py --n-values 50 100 150 200 250 300 --runs 5
```

Esto genera:
- Archivos de simulación en `data/sim_<N>N_<timestamp>.txt`
- Archivo de tiempos en `data/timing.txt`

---

## Parámetros Físicos

| Parámetro            | Valor     | Unidad |
|----------------------|-----------|--------|
| Diámetro recinto (L) | 80        | m      |
| Radio recinto        | 40        | m      |
| Radio obstáculo (r₀) | 1         | m      |
| Radio partícula (r)  | 1         | m      |
| Masa partícula (m)   | 1         | kg     |
| Velocidad inicial (v₀)| 1        | m/s    |
| Tiempo final (t_f)   | 5         | s      |

---

## Generación de Gráficos (Python)

### Modo archivo único (Incisos 1.3 y 1.4)

```bash
python graphics/plot_metrics.py data/sim_100N_20260410_1530.txt
```

### Modo directorio completo (Incisos 1.1–1.4)

```bash
python graphics/plot_metrics.py data/
```

### Scripts individuales

Cada inciso tiene su script independiente:

---

### Inciso 1.1 — Tiempo de ejecución vs N

```bash
python graphics/plot_execution_time.py
```

**Requisito:** Archivo `data/timing.txt` generado por `run_batch.py`.

**Genera:** `graphics/output/inciso_1_1_execution_time.png`

---

### Inciso 1.2 — Scanning rate ⟨J⟩(N)

```bash
python graphics/plot_scanning_rate.py
```

Lee todos los archivos `sim_*.txt` de `data/`, agrupa por N, calcula C_fc(t) y obtiene J como la pendiente de la interpolación lineal.

**Genera:**
- `graphics/output/inciso_1_2_scanning_rate.png` — ⟨J⟩ vs N con barras de error
- `graphics/output/inciso_1_2_cfc_curve.png` — C_fc(t) para el mayor N

---

### Inciso 1.3 — Fracción de partículas usadas F_u(t)

```bash
python graphics/plot_fraction_used.py data/sim_100N_20260410_1530.txt
```

Muestra la evolución temporal y el régimen estacionario.

**Genera:** `graphics/output/inciso_1_3_fraction_used_N<N>.png`

---

### Inciso 1.4 — Perfiles radiales

```bash
python graphics/plot_radial_profiles.py data/sim_100N_20260410_1530.txt
```

Calcula perfiles de densidad, velocidad normalizada y flujo para partículas frescas moviéndose hacia el centro (R_j · v_j < 0), en capas concéntricas de dS = 0.2 m.

**Genera:** `graphics/output/inciso_1_4_radial_profiles_N<N>.png`

---

### Animación

```bash
python graphics/animate_system1.py data/sim_100N_20260410_1530.txt
python graphics/animate_system1.py data/sim_100N_20260410_1530.txt --save anim.gif --skip 5
```

---

## Formato del Archivo de Salida

```
# N=100 L=80.0 R_enclosure=40.0 r0=1.0 r=1.0 m=1.0 v0=1.0 t_final=5.0
# FORMAT: SNAPSHOT lines start with 'S', followed by time,
# then N lines of: id x y vx vy state(F/U)
# EVENT lines start with 'E': time type id1 [id2]
S 0.000000e+00
0 1.234567e+01 -5.678901e+00 7.071068e-01 7.071068e-01 F
1 -2.345678e+01 1.234567e+01 -5.000000e-01 8.660254e-01 U
...
S 1.000000e-02
...
```

- **S**: línea de snapshot con el tiempo del evento
- Cada partícula: `id x y vx vy estado`
- Estado: `F` (fresca/verde) o `U` (usada/violeta)
- Unidades MKS, notación científica

---

## Flujo de Trabajo Completo

```bash
# 1. Compilar el motor Java
cd engine && mvn compile && cd ..

# 2. Ejecutar batch de simulaciones
python graphics/run_batch.py --n-values 50 100 150 200 250 300 --runs 5

# 3. Generar todos los gráficos
python graphics/plot_metrics.py data/

# 4. (Opcional) Animación
python graphics/animate_system1.py data/sim_300N_*.txt --skip 10
```

---

## Arquitectura

```
┌─────────────────────────┐         ┌─────────────────────────┐
│      Motor Java         │         │   Post-procesamiento    │
│  (EventDrivenSimulation)│  ─────> │       Python            │
│                         │ archivo │                         │
│  • PriorityQueue<Event> │  .txt   │  • plot_metrics.py      │
│  • Colisión-colisión    │         │  • matplotlib + numpy   │
│  • Colisión-pared       │         │                         │
│  • Colisión-obstáculo   │         │  Gráficos:              │
│  • Invalidación por     │         │  • 1.1: T(N)            │
│    conteo de colisiones │         │  • 1.2: <J>(N)          │
│                         │         │  • 1.3: F_u(t)          │
│  Output: data/*.txt     │         │  • 1.4: ρ(S), v(S),    │
│                         │         │         J_in(S)         │
└─────────────────────────┘         └─────────────────────────┘
```

---

## Matemáticas de Colisión

### Tiempo de colisión partícula-partícula (Slide 14/20)
- Δr = r_j - r_i, Δv = v_j - v_i
- σ = r_i + r_j
- d = (Δv·Δr)² - (Δv·Δv)(Δr·Δr - σ²)
- t_c = -(Δv·Δr + √d) / (Δv·Δv) si Δv·Δr < 0 y d ≥ 0

### Post-choque elástico (Slide 36/37)
- J = 2·m_i·m_j·(Δv·Δr) / (σ·(m_i + m_j))
- v_i' = v_i + (J/m_i)·n̂, v_j' = v_j - (J/m_j)·n̂

### Colisión con paredes circulares
- |r + v·t|² = R_eff² → ecuación cuadrática en t
- Pared exterior: R_eff = R_enclosure - r → raíz positiva mayor
- Obstáculo: R_eff = r_0 + r → raíz positiva menor, si r·v < 0
