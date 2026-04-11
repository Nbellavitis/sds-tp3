# TP3 — Dinámica Molecular Dirigida por Eventos

## Sistema 1: Scanning Rate en Recinto Cerrado con Obstáculo Fijo

Simulación de partículas en un recinto circular con un obstáculo fijo en el centro, utilizando Dinámica Molecular Dirigida por Eventos (Event-Driven Molecular Dynamics).

---

## Estructura del Proyecto

```
sds-tp3/
├── engine/                          # Motor de simulación (Java 17+)
│   ├── pom.xml
│   └── src/main/java/ar/edu/itba/sds/tp3/
│       ├── Particle.java            # Partícula con colisiones y estados
│       ├── Event.java               # Eventos para la PriorityQueue
│       └── EventDrivenSimulation.java  # Motor principal (main)
│
├── graphics/                        # Post-procesamiento (Python 3)
│   ├── plot_metrics.py              # Script principal (orquestador de todos los incisos)
│   ├── plot_execution_time.py       # Inciso 1.1: Tiempo de ejecución vs N
│   ├── plot_scanning_rate.py        # Inciso 1.2: Scanning rate <J>(N) + C_fc(t)
│   ├── plot_fraction_used.py        # Inciso 1.3: F_u(t)
│   ├── plot_radial_profiles.py      # Inciso 1.4: ρ(S), v(S), J_in(S) + vs N en S≈2
│   ├── animate_system1_mru.py       # Animación interpolada con MRU entre eventos
│   └── run_batch.py                 # Runner de batch para múltiples N y seeds
│
├── data/                            # Archivos de salida de la simulación
│   └── sim_<N>N_<timestamp>_s<seed>.txt
│   └── cache/                       # Cachés de análisis generadas por Java para post-procesado rápido
│
└── README.md
```

---

## Requisitos

### Java
- **JDK 17+** (probado con JDK 21)

### Python
- **Python 3.8+**
- **Dependencias:** `pip install numpy matplotlib`

---

## Compilación (Java)

```bash
cd engine
# Opción 1: Con Maven (si está instalado)
mvn compile

# Opción 2: Con javac directamente
cd ..
javac -d engine/target/classes \
  engine/src/main/java/ar/edu/itba/sds/tp3/Particle.java \
  engine/src/main/java/ar/edu/itba/sds/tp3/Event.java \
  engine/src/main/java/ar/edu/itba/sds/tp3/EventDrivenSimulation.java
```

---

## Ejecución de la Simulación

### Ejecución simple

```bash
java -cp engine/target/classes ar.edu.itba.sds.tp3.EventDrivenSimulation -N 200
java -cp engine/target/classes ar.edu.itba.sds.tp3.EventDrivenSimulation -N 200 -t_final 8.0
```

### Parámetros

| Parámetro | Descripción                     | Default |
|-----------|---------------------------------|---------|
| `-N`      | Número de partículas            | 100     |
| `-seed`   | Semilla para reproducibilidad   | aleatorio |
| `-runs`   | Número de realizaciones         | 1       |
| `-t_final` | Tiempo final de simulación [s] | 5.0     |

### Múltiples realizaciones

```bash
java -cp engine/target/classes ar.edu.itba.sds.tp3.EventDrivenSimulation -N 200 -runs 5 -seed 42 -t_final 8.0
```

### Batch (múltiples N)

```bash
python graphics/run_batch.py --n-values 50 100 150 200 250 300 --runs 5
python graphics/run_batch.py --n-start 50 --n-stop 500 --n-step 50 --runs 5
python graphics/run_batch.py --for-1-1 --runs 5
python graphics/run_batch.py --for-1-1 --n-max 650 --runs 5
python graphics/run_batch.py --n-values 50 100 150 200 250 300 --runs 5 --t-final 8.0
```

---

## Parámetros Físicos

| Parámetro            | Valor | Unidad |
|----------------------|-------|--------|
| Diámetro recinto (L) | 80    | m      |
| Radio recinto        | 40    | m      |
| Radio obstáculo (r₀) | 1     | m      |
| Radio partícula (r)  | 1     | m      |
| Masa partícula (m)   | 1     | kg     |
| Velocidad inicial (v₀)| 1    | m/s    |
| Tiempo final (t_f)   | 5 (default, configurable) | s |

---

## Generación de Gráficos

Los scripts de `graphics/` ya no parsean directamente los `sim_*.txt` grandes en cada corrida.
En la primera ejecución construyen una caché compacta en `data/cache/` usando Java, y luego Python solo carga esos datos y grafica.

### Todos los incisos (desde un directorio con múltiples archivos)

```bash
python graphics/plot_metrics.py data/
```

### Incisos 1.3 y 1.4 (desde un archivo individual)

```bash
python graphics/plot_metrics.py data/sim_300N_20260409_235938_s42.txt
```

---

## Detalle por Inciso

### Inciso 1.1 — Tiempo de ejecución vs N

```bash
python graphics/run_batch.py --for-1-1 --runs 5
python graphics/run_batch.py --for-1-1 --n-max 650 --runs 5
python graphics/plot_execution_time.py
```

- **Requisito:** `data/timing_1_1.txt` (generado por `run_batch.py --for-1-1`)
- La consigna de 1.1 usa `t_f = 5 s`
- `--for-1-1` fija automáticamente `N = 50, 100, ..., Nmax` con paso `50` y `t_f = 5 s`
- `--n-max` te deja elegir ese `N` máximo; por default vale `500`
- El gráfico de 1.1 solo usa archivos de timing con `t_f = 5 s` (si no, avisa y no grafica)
- **Genera:** `graphics/output/inciso_1_1_execution_time.png`

---

### Inciso 1.2 — Scanning rate ⟨J⟩(N)

```bash
python graphics/plot_scanning_rate.py
```

- Lee las líneas `E <time> <id>` de los archivos de simulación (transiciones F→U exactas)
- Construye C_fc(t) y calcula J como la pendiente de la interpolación lineal
- **Genera:**
  - `graphics/output/inciso_1_2_scanning_rate.png` — ⟨J⟩ vs N con error
  - `graphics/output/inciso_1_2_cfc_curve.png` — C_fc(t) con ajuste lineal

---

### Inciso 1.3 — Fracción de partículas usadas F_u(t)

```bash
python graphics/plot_fraction_used.py data/sim_300N_20260409_235938_s42.txt   # un archivo
python graphics/plot_fraction_used.py data/                  # directorio
```

- Muestra `F_u(t)` para todas las realizaciones de cada `N`, superpuestas en un mismo gráfico
- Marca con líneas punteadas los valores manuales de `t_est` y `F_est` para cada `N`
- **Genera:**
  - `graphics/output/inciso_1_3_fraction_used_N<N>.png` — `F_u(t)` para cada `N`, con realizaciones superpuestas
  - `graphics/output/inciso_1_3_t_est_vs_N.png` — valores manuales de `t_est` en función de `N`
  - `graphics/output/inciso_1_3_fest_vs_N.png` — valores manuales de `F_est` en función de `N`

---

### Inciso 1.4 — Perfiles radiales

```bash
python graphics/plot_radial_profiles.py data/sim_300N_20260409_235938_s42.txt  # un archivo
python graphics/plot_radial_profiles.py data/                 # directorio
```

- S = distancia radial desde el centro (origen)
- Capas concéntricas dS = 0.2 m
- Solo partículas frescas con R·v < 0 (velocidad apuntando al centro)
- v_f^in = (R·v)/|R| (velocidad radial proyectada)
- **Genera:**
  - `graphics/output/inciso_1_4_radial_profiles_N<N>.png` — ρ(S), |v(S)|, J_in(S)
  - `graphics/output/inciso_1_4_S2_vs_N.png` — ρ, v, J_in en S≈2 vs N

---

### Animación

```bash
python graphics/animate_system1_mru.py data/sim_300N_20260409_235938_s42.txt
python graphics/animate_system1_mru.py data/sim_300N_20260409_235938_s42.txt --save anim.gif --fps 20 --speed 2.0
```

---

## Formato del Archivo de Salida

```
# N=300 L=80.0 R_enclosure=40.0 r0=1.0 r=1.0 m=1.0 v0=1.0 t_final=5.0
# FORMAT: SNAPSHOT lines start with 'S', followed by time,
# then N lines of: id x y vx vy state(F/U)
# EVENT lines start with 'E': time id
S 0.000000e+00
0 1.234567e+01 -5.678901e+00 7.071068e-01 7.071068e-01 F
1 -2.345678e+01 1.234567e+01 -5.000000e-01 8.660254e-01 F
...
E 2.619697e+00 196       <-- F->U transition (for C_fc counting)
S 2.619697e+00
...
```

- **`S <time>`**: snapshot en tiempos de evento — seguido de N líneas con estado de cada partícula
- **`E <time> <id>`**: evento F→U — partícula fresca contactó el obstáculo central
- Estado: `F` (fresca/verde) o `U` (usada/violeta)
- Unidades MKS, notación científica

---

## Flujo de Trabajo Completo

```bash
# 1. Compilar
javac -d engine/target/classes engine/src/main/java/ar/edu/itba/sds/tp3/*.java

# 2. Batch de simulaciones
python graphics/run_batch.py --n-values 50 100 150 200 250 300 --runs 5
# opcional: cambiar el tiempo final de simulación
python graphics/run_batch.py --n-values 50 100 150 200 250 300 --runs 5 --t-final 8.0

# 3. Todos los gráficos
python graphics/plot_metrics.py data/

# 4. Animación (opcional)
python graphics/animate_system1_mru.py data/sim_300N_20260409_235938_s42.txt --fps 20 --speed 2.0
```

---

## Matemáticas de Colisión (basadas en Teórica 3)

### Tiempo de colisión partícula-partícula (Diapositiva 14)
```
Δr = r_j - r_i,  Δv = v_j - v_i,  σ = R_i + R_j
d = (Δv·Δr)² - (Δv·Δv)(Δr·Δr - σ²)
t_c = -(Δv·Δr + √d) / (Δv·Δv)   si Δv·Δr < 0 y d ≥ 0
```

### Post-choque elástico (Diapositiva 20)
```
J = 2·m_i·m_j·(Δv·Δr) / (σ·(m_i + m_j))
Jx = J·Δx/σ,  Jy = J·Δy/σ
v_i' = v_i + (J/m_i)·n̂,  v_j' = v_j - (J/m_j)·n̂
```

### Colisión con paredes circulares
```
|r + v·t|² = R_eff²  →  ecuación cuadrática en t
Pared exterior: R_eff = R_enclosure - r  →  raíz positiva (la menor > 0)
Obstáculo:      R_eff = r₀ + r          →  raíz positiva menor
```

### Colisión obstáculo fijo (Diapositivas 35-37)
```
Reflexión especular: v' = v - 2(v·n̂)n̂
n̂ = r/|r| (normal saliente del obstáculo)
```
