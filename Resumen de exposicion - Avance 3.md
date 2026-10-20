---
tags:
  - metodos-numericos
  - CC2104
  - exposicion
  - avance3
titulo: Resumen de exposición — Avance 3
duracion_estimada: 5-7 min
vinculado:
  - "[[5.3 Avance 3 - Implementación y Borrador de Resultados]]"
  - "[[MN P2 - Esquema Grupo 4]]"
  - "[[MN P2 - Esquema.canvas]]"
---

# 🎤 Resumen para exponer el Avance 3

> [!info] Objetivo de la exposición
> Dejar claro que este avance = **"corrimos el código + mostramos resultados preliminares + empezamos a discutir"**, y que falta el **análisis paramétrico completo** y el **error vs. referencia** para la entrega final.

---

## ⏱ Guion rápido (5–7 min)

### 0:00 · Apertura (una frase)
> *"El Avance 3 es el paso previo a la entrega final: implementamos los métodos numéricos en MATLAB, corrimos el caso base y el barrido de β, y ya estamos esbozando la discusión."*

### 0:30 · El problema
- Ecuación de **advección–difusión 1D transitoria** con **difusividad térmica variable** α(T).
- Flujo sobre un **tubo cilíndrico**: desde la corriente libre (x=0) hasta la pared (x=L).
- **β ≠ 0 hace la ecuación no lineal** → no hay solución exacta → métodos numéricos.

### 1:15 · La metodología (la cadena de 4 métodos)
> *"Usamos cuatro métodos, uno dentro de otro."*

1. **Diferencias finitas** → discretiza el **espacio** (difusión conservativa en puntos medios)
2. **SDIRK4** (5 etapas) → avanza el **tiempo**; implícito ⇒ permite h mayores
3. **Newton–Raphson** → resuelve lo **no lineal**; el Jacobiano es **tridiagonal**
4. **Gauss–Seidel** → resuelve el **sistema lineal** dentro de cada iteración de Newton

**Criterios de fallo:** Newton máx. 20 iter · Gauss–Seidel máx. 500 iter.

### 2:15 · Caso base
| Parámetro | Valor |
|---|---|
| L · N · Δx | 1 m · 51 nodos · 0,02 m |
| v · α₀ · β | 0,5 m/s · 0,01 m²/s · 0,2 |
| h · pasos | 0,01 s × 50 → **t = 0,5 s** |

**Frente x_f** = punto donde T = 50 °C (valor medio), por interpolación lineal.

### 3:00 · Resultados
**Caso base**
- Frente **x_f = 0,4739 m**
- Newton: **787** iteraciones (≈3,15 por etapa)
- Gauss–Seidel: **9 140** (≈11,6 por cada Newton)
- **0 etapas sin convergencia** · sin oscilaciones · el frente **no llega** a la pared

**Barrido de β**
| β | x_f (m) | G-S | Tiempo (s) |
|---|---|---|---|
| 0 | 0,2670 | 1 925 | 0,70 |
| 0,1 | 0,3965 | 6 372 | 2,13 |
| 0,2 | 0,4739 | 9 140 | 3,39 |
| 0,4 | 0,5882 | 14 294 | 5,07 |

### 4:00 · Los 3 hallazgos (lo importante de decir en voz alta)

> [!quote] 1 — El esquema es robusto
> **Cero etapas sin convergencia** con cualquier β. El método funciona.

> [!quote] 2 — Newton apenas depende de la no linealidad
> De β=0 a β=0,1 sube de 500 a 774 iteraciones, pero de 0,1 a 0,4 casi no cambia (774 → 811). **El costo creciente NO lo causa Newton.**

> [!quote] 3 — El costo lo causa Gauss–Seidel
> Iteraciones de G-S por cada Newton: **3,9 → 8,2 → 11,6 → 17,6**. El tiempo se multiplica **×7,2** (0,70 → 5,07 s). *Ahí está el cuello de botella.*

**Extra (física):** subir β aleja el frente **+120 %** respecto al caso lineal, porque α crece con T y la zona caliente conduce mejor. Y con β=0, x_f ≈ 0,267 m ≈ la advección pura `v·t = 0,25 m`.

### 4:45 · Discusión vs. objetivos
| Objetivo | Estado |
|---|---|
| **O1** Implementar DF+SDIRK4+NR+GS | ✅ hecho, 0 fallas |
| **O2** Error vs. solución de referencia | ⏳ **pendiente** |
| **O3** Robustez: β, v, Δx, h | 🟡 **solo β**; faltan v, Δx, h |
| **O4** Posición del frente | 🟡 hecho, pero **falta obtenerlo por código** en cada instante |

### 5:15 · Qué falta para la entrega final
- [ ] **Error máximo** vs. referencia (O2)
- [ ] Barridos de **v, Δx, h** (O3)
- [ ] **x_f por código** en cada instante, no leyéndolo de la figura (O4)
- [ ] Formato **APA** en figuras y tablas
- [ ] Detalle **computacional** de la metodología en el paper

### 5:45 · Cierre
> *"En síntesis: el esquema convergió siempre, el frente avanza más lejos a mayor β, y el costo de cómputo lo domina Gauss–Seidel. Lo que queda es completar el análisis paramétrico y cuantificar el error."*

---

## 🗣 Frases clave (para no trabarse)

- *"La no linealidad entra por β; si β=0 recuperamos el caso lineal."*
- *"El Jacobiano es tridiagonal porque cada nodo solo depende de sus vecinos."*
- *"SDIRK4 es implícito, por eso tolera pasos de tiempo mayores."*
- *"Cada etapa de Runge-Kutta exige resolver una ecuación no lineal, y ahí entran Newton-Raphson y Gauss-Seidel."*
- *"El cuello de botella es el solver lineal, no Newton."*

---

## ❓ Preguntas probables y respuestas

**¿Por qué SDIRK4 y no RK4 clásico?**
> Porque RK4 es **explícito** y no necesita resolver sistemas; SDIRK4 es **implícito** y permite **h mayores** y es estable para cambios bruscos de temperatura (Wang 2015; Deng 2014; Levý & May 2022). Además no lo vimos en el curso, pero su orden 4 lo justifica.

**¿Por qué Newton-Raphson converge tan rápido?**
> Porque la aproximación inicial está cerca de la solución (Deuflhard, 1974) y el Jacobiano es **tridiagonal**, así que cada iteración es barata.

**¿Qué mide β?**
> Cuánto cambia la difusividad con la temperatura: `α(T) = α₀[1 + β(T − T₀)]`. Con β=0,2 y T=100 °C, α = 21 α₀.

**¿Cómo se define el frente?**
> Punto donde T alcanza el valor medio **(T₁+T_N)/2 = 50 °C**: primer nodo con T ≤ 50 + interpolación lineal con el anterior.

**¿Por qué el tiempo crece tanto con β?**
> No es Newton (casi constante), son las **iteraciones de Gauss–Seidel** que suben de 3,9 a 17,6 por cada Newton.

**¿Qué falta?**
> El error vs. la solución de referencia, y los barridos de v, Δx y h.

---

## 📋 Diapositivas sugeridas

1. **Título + problema** (ecuación y por qué es no lineal)
2. **Cadena de métodos** (diagrama DF → SDIRK4 → NR → GS)
3. **Caso base** (Tabla 1)
4. **Resultados caso base** (Tabla 2 + Figura 1)
5. **Evolución del perfil** (Figura 2 + x_f en t=0,1…0,5)
6. **Barrido de β** (Tabla 3)
7. **Discusión ↔ objetivos** (estado O1–O4)
8. **Pendientes + cierre**

> [!tip] Consejo
> Si solo tienes **3 minutos**: problema (30 s) → cadena de métodos (45 s) → Tabla 3 + *"el costo lo domina Gauss–Seidel"* (60 s) → pendientes (45 s).

---

#CC2104 #avance3 #exposicion #guion
