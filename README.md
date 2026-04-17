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
│   ├── run_batch.py                 # Runner de batch para múltiples N y seeds
│   └── run_animation_sim.py         # Runner separado para generar salidas densas para animación
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
| `-snapshot_every_events` | Guarda una snapshot cada K colisiones (para animación usar K=1 o K bajo) | 10 |
| `--no-output` | No escribe `sim_*.txt`; útil para medir tiempo sin costo de I/O | desactivado |

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
- La consigna de 1.1 usa los primeros `5 s` simulados
- `--for-1-1` sigue siendo una forma cómoda de generar el batch dedicado, pero ya no hace falta que la simulación completa termine exactamente en `t_f = 5 s`
- Ese batch corre la simulación con `--no-output`: genera solo `data/timing_1_1.txt`, no crea `sim_*.txt`, y así la medición no arrastra costo de I/O
- `data/timing.txt` guarda tiempo total de realización; `data/timing_1_1.txt` guarda solo el tiempo hasta los primeros `5 s` simulados
- `--n-max` te deja elegir ese `N` máximo; por default vale `500`
- El gráfico de 1.1 acepta archivos de timing medidos hasta `t=5 s` simulados, aun si la corrida completa tuvo `t_f > 5 s`
- Cada punto grafica el tiempo medio de ejecución hasta esos primeros `5 s` simulados para ese `N`
- Las barras verticales representan la desviación estándar poblacional entre realizaciones:
  promedio: `T_prom(N) = (1/R) * sum_r T_r(N)`
  error: `sigma_T(N) = sqrt((1/R) * sum_r (T_r(N) - T_prom(N))^2)`
- El gráfico no incluye ajuste exponencial
- **Genera:** `graphics/output/inciso_1_1_execution_time.png`

---

### Inciso 1.2 — Scanning rate ⟨J⟩(N)

```bash
python graphics/plot_scanning_rate.py
```

- Lee las líneas `E <time> <id>` de los archivos de simulación (transiciones F→U exactas)
- Construye C_fc(t) y calcula J como la pendiente de la interpolación lineal
- Para cada `N`, el punto mostrado es el promedio de los `J_r` de todas las realizaciones
- Las barras verticales representan la desviación estándar poblacional entre realizaciones:
  promedio: `J_prom(N) = (1/R) * sum_r J_r(N)`
  error: `sigma_J(N) = sqrt((1/R) * sum_r (J_r(N) - J_prom(N))^2)`
- El gráfico `C_fc(t)` se hace para el mayor `N` cargado y aclara dentro de la figura que las líneas sólidas son los contactos acumulados y las punteadas el ajuste lineal
- **Genera:**
  - `graphics/output/inciso_1_2_scanning_rate.png` — ⟨J⟩ vs N con barras verticales de dispersión entre realizaciones
  - `graphics/output/inciso_1_2_cfc_curve.png` — C_fc(t) con ajuste lineal

---

### Inciso 1.3 — Fracción de partículas usadas F_u(t)

```bash
python graphics/plot_fraction_used.py data/sim_300N_20260409_235938_s42.txt   # un archivo
python graphics/plot_fraction_used.py data/                  # directorio
```

- En modo directorio, genera las figuras temporales solo para `N = 100, 300, 500, 800` (si están disponibles)
- Muestra `F_u(t)` para todas las realizaciones de cada `N`, superpuestas en un mismo gráfico
- Marca con línea punteada vertical el `t_est` manual de cada `N`
- Calcula `F_est` promediando los valores de `F_u(t)` a partir de `t_est`
- La leyenda distingue realizaciones, pero ya no muestra la `seed`
- En `1.3` no se grafican barras de error: `t_est` sale de la tabla manual `MANUAL_T_EST_BY_N` y `F_est` se calcula desde los datos cargados
- **Genera:**
  - `graphics/output/inciso_1_3_fraction_used_N<N>.png` — `F_u(t)` para cada `N` seleccionado, con realizaciones superpuestas
  - `graphics/output/inciso_1_3_t_est_vs_N.png` — valores manuales de `t_est` en función de `N`
  - `graphics/output/inciso_1_3_fest_vs_N.png` — valores calculados de `F_est` en función de `N`

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
- En modo directorio, genera los perfiles radiales solo para `N = 100, 300, 500, 800` (si están disponibles)
- Los perfiles se calculan usando únicamente las snapshots guardadas
- Para cada corona, la velocidad radial media se calcula acumulando todas las partículas observadas en esa corona en todas las snapshots y realizaciones, y recién ahí se promedia
- `ρ(S)` y `J_in(S)` se construyen a partir de esos acumuladores por snapshot; `J_in(S) = ρ(S) * |v(S)|`
- Las barras verticales siguen mostrando la desviación estándar poblacional entre realizaciones
- En los gráficos de `S≈2`, se toma la capa cuyo centro es el primero que cumple `S_c ≥ 2 m`
- **Genera:**
  - `graphics/output/inciso_1_4_rho_profile_N<N>.png` — perfil radial de `ρ(S)` para cada `N` seleccionado
  - `graphics/output/inciso_1_4_velocity_profile_N<N>.png` — perfil radial de `|v(S)|` para cada `N` seleccionado
  - `graphics/output/inciso_1_4_flux_profile_N<N>.png` — perfil radial de `J_in(S)` para cada `N` seleccionado
  - `graphics/output/inciso_1_4_S2_rho_vs_N.png` — `ρ` en `S≈2` vs `N`
  - `graphics/output/inciso_1_4_S2_velocity_vs_N.png` — `|v|` en `S≈2` vs `N`
  - `graphics/output/inciso_1_4_S2_flux_vs_N.png` — `J_in` en `S≈2` vs `N`

---

### Animación

```bash
python graphics/animate_system1_mru.py data/sim_300N_20260409_235938_s42.txt
python graphics/animate_system1_mru.py data/sim_300N_20260409_235938_s42.txt --save anim.gif --fps 20 --speed 2.0

# runner separado para generar archivos orientados a animación
python graphics/run_animation_sim.py --n 300 --seed 42 --t-final 8.0 --snapshot-every-events 1
```

---

## Formato del Archivo de Salida

```
# N=300 L=80.0 R_enclosure=40.0 r0=1.0 r=1.0 m=1.0 v0=1.0 t_final=5.0 snapshot_every_events=1
# FORMAT: SNAPSHOT lines start with 'S', followed by time,
# then N lines of: id x y vx vy state(F/U)
# EVENT lines start with 'E': time id
S 0.000000e+00
0 1.234567e+01 -5.678901e+00 7.071068e-01 7.071068e-01 F
1 -2.345678e+01 1.234567e+01 -5.000000e-01 8.660254e-01 F
...
E 2.619697e+00 196       <-- F->U transition (for C_fc counting)
...
S 2.619697e+00          <-- snapshot guardada cada K eventos (K configurable), mas t=0 y el estado final
...
```

- **`S <time>`**: snapshot guardada en `t=0`, cada `K` colisiones procesadas y al final (con `K = snapshot_every_events`) — seguida de N líneas con estado de cada partícula
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
