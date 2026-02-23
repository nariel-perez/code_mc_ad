# Informe de bugs y correcciones — funciona_ok

Este informe resume los bugs detectados en `funciona_ok` (según la revisión de 2026) y cómo se corrigieron. Se mantuvo **sin cambios** el array `LOCATE(5000,10)` según lo solicitado.

---

## 1. Bugs encontrados y correcciones aplicadas

### 1.1 Estructura.f90 — Compresión incompleta de `EPSAC` y `SGC`

- **Problema**: al filtrar átomos dentro de la celda, `EPSAC` y `SGC` quedaban con índices no comprimidos.
- **Riesgo**: potencial superficie–adsorbato calculado con parámetros incorrectos.
- **Corrección**: se copian `EPSAC(i)` y `SGC(i)` a `EPSAC(imax)`/`SGC(imax)` junto a `RXC/RYC/RZC/QAC`.

### 1.2 Potencial.f90 / PotencialFF.f90 — Lectura inconsistente de `LJ.dat`

- **Problema**: `Potencial` leía 1 valor en la primera línea y `PotencialFF` leía 2.
- **Riesgo**: fallo de lectura en un formato u otro y tablas `USS` en cero.
- **Corrección**: lectura robusta de la primera línea con buffer, aceptando 1 o 2 valores. Se define un valor por defecto para `REDELEC` si no viene en el archivo.

### 1.3 Out.f90 — Eliminación sin protección cuando `N(MOLKIND)=0`

- **Problema**: posible acceso fuera de rango de `LOCATE` al eliminar moléculas inexistentes.
- **Riesgo**: segfault o lectura de memoria basura.
- **Corrección**: guard clause inmediata si `N(MOLKIND) <= 0`.

### 1.4 Main.f90 — `Z(NMOLEC)` sobrescrito dentro del bucle

- **Problema**: la actividad de la última especie se sobrescribía en cada iteración.
- **Riesgo**: actividad fija (55.55) aunque no corresponda a la especie.
- **Corrección**: se movió la asignación fuera del bucle y se dejó comentario explícito.

### 1.5 Potin.f90 / Potout.f90 — Uso de `WIJ` no inicializada

- **Problema**: `WIJ` se sumaba sin asignación previa.
- **Riesgo**: acumulación de basura si se usa virial en el futuro.
- **Corrección**: se eliminó la suma de `WIJ` y variables asociadas.

### 1.6 Potencial.f90 — `DELTW` no inicializada por punto de grilla

- **Problema**: `DELTW` se acumulaba sin reinicio por cada punto.
- **Riesgo**: valores residuales si se habilita virial.
- **Corrección**: inicialización de `DELTW = 0.0` dentro del bucle interno.

### 1.7 Ranf.f90 — Generador con período muy corto

- **Problema**: LCG con período ~1e6.
- **Riesgo**: correlaciones y repetición de la secuencia.
- **Corrección**: reemplazo por `RANDOM_NUMBER` manteniendo la misma interfaz `RANF(DUMMY)`.

### 1.8 RotationModule.f90 — PI de baja precisión y sin clamp de índices

- **Problema**: `PI` truncado y riesgo de índice fuera de rango.
- **Riesgo**: acceso inválido a `TABLE_SEN/TABLE_COS`.
- **Corrección**: `PI = ACOS(-1.0)` y clamp `MIN/MAX` para índices.

### 1.9 Main.f90 — `Ener###.TXT` abierto/cerrado sin uso

- **Problema**: apertura/cierre por iteración sin escritura.
- **Riesgo**: I/O innecesario.
- **Corrección**: se eliminó el `open/close` del unit 51.

---

## 2. Cambios NO aplicados por solicitud

- `LOCATE(5000,10)` permanece igual en `SimulationData.f90`.

---

## 3. Archivos modificados

- `Estructura.f90`
- `Potencial.f90`
- `PotencialFF.f90`
- `Out.f90`
- `Main.f90`
- `Potin.f90`
- `Potout.f90`
- `RotationModule.f90`
- `Ranf.f90`

