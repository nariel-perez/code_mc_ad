# Revisión de bugs y problemas — modulos_f90 (febrero 2026)

Análisis exhaustivo del código GCMC tras lectura completa de todos los archivos fuente.
Se clasifica cada hallazgo con severidad y se incluyen ejemplos concretos de cómo falla.

---

## PARTE 1: BUGS Y PROBLEMAS DETECTADOS

---

### BUG 1 — CRÍTICO: `Estructura.f90` no comprime `EPSAC` ni `SGC` al filtrar átomos

**Archivo:** `Estructura.f90`, líneas 65–88

**Qué hace el código:**
Lee todos los átomos de la superficie desde el archivo `nam`. Para cada átomo, verifica si está dentro de la celda de simulación (coordenadas fraccionales en [-0.5, 0.5]). Los que pasan el filtro se copian a los arrays `RXC`, `RYC`, `RZC`, `QAC` con un índice comprimido `imax`.

**El problema:**
Los arrays `EPSAC` y `SGC` (epsilon y sigma de cada átomo de la superficie) **no se comprimen**. Se leen en el bucle (línea 67) con índice `i` (1..NC original), pero después del filtrado, `Potencial.f90` los recorre con índice `J = 1..NC` (NC ya actualizado al valor filtrado). Los valores en `EPSAC(J)` y `SGC(J)` corresponden a los primeros J átomos del archivo original, no a los átomos que pasaron el filtro.

**Ejemplo concreto:**

Supongamos que el archivo de superficie tiene 5 átomos:

```
Átomo 1: x=1.0  y=2.0  z=3.0   eps=28.0  sigma=3.4   (DENTRO de la celda)
Átomo 2: x=50.0 y=50.0 z=50.0  eps=80.0  sigma=3.05  (FUERA de la celda)
Átomo 3: x=2.0  y=1.0  z=1.0   eps=36.0  sigma=3.31  (DENTRO de la celda)
Átomo 4: x=60.0 y=60.0 z=60.0  eps=137.0 sigma=3.59  (FUERA de la celda)
Átomo 5: x=3.0  y=3.0  z=2.0   eps=202.0 sigma=3.83  (DENTRO de la celda)
```

Después del filtrado:
- `imax = 3` (solo 3 átomos dentro)
- `RXC(1), RYC(1), RZC(1)` = coords del átomo 1 ✓
- `RXC(2), RYC(2), RZC(2)` = coords del átomo 3 ✓
- `RXC(3), RYC(3), RZC(3)` = coords del átomo 5 ✓
- `QAC(1)` = carga del átomo 1 ✓
- `QAC(2)` = carga del átomo 3 ✓
- `QAC(3)` = carga del átomo 5 ✓

Pero:
- `EPSAC(1) = 28.0` (átomo 1) ✓
- `EPSAC(2) = 80.0` (átomo 2, ¡que está FUERA!) ✗ debería ser 36.0 (átomo 3)
- `EPSAC(3) = 36.0` (átomo 3) ✗ debería ser 202.0 (átomo 5)
- `SGC(1) = 3.4` (átomo 1) ✓
- `SGC(2) = 3.05` (átomo 2, ¡FUERA!) ✗ debería ser 3.31
- `SGC(3) = 3.31` (átomo 3) ✗ debería ser 3.83

**Consecuencia:** `Potencial.f90` calcula las interacciones adsorbato-superficie con parámetros de Lennard-Jones incorrectos. El mapa de energía `UADS` está mal calculado, lo que afecta toda la simulación: aceptación/rechazo de inserciones, cálculo de isotermas, calores isostéricos, etc.

**Dónde falla exactamente — Potencial.f90, líneas 73-74:**

```fortran
SIGMA1 = SIGMA * (SGC(J) + SIGM(IPOT)) / (2 * SIGMETANO)   ! SGC(J) es incorrecto
FACTOR = SQRT(EPSI(IPOT) * EPSAC(J)) / EPS                  ! EPSAC(J) es incorrecto
```

**Corrección necesaria:** Dentro del bloque `IF (all(abs(s) <= 0.5_rk)) THEN`, agregar:

```fortran
EPSAC(imax) = EPSAC(i)
SGC(imax)   = SGC(i)
```

de forma análoga a como ya se hace con `QAC(imax) = QAC(i)`.

---

### BUG 2 — CRÍTICO: Lectura inconsistente de `LJ.dat` entre `PotencialFF.f90` y `Potencial.f90`

**Archivos:** `PotencialFF.f90` línea 31, `Potencial.f90` línea 35

**Qué hace el código:**
Ambas subrutinas abren y leen el archivo `LJ.dat`. La primera línea del archivo contiene información sobre los tipos de átomos.

**El problema:**
`PotencialFF.f90` lee **dos valores** de la primera línea:

```fortran
READ(11, *) NKIND, REDELEC     ! espera: "3  0.35" (por ejemplo)
```

`Potencial.f90` lee **un solo valor** de la primera línea:

```fortran
READ(11, *) NKIND              ! espera: "3"
```

Si el archivo `LJ.dat` tiene el formato que espera `PotencialFF`:

```
3  0.35
1  28.0  3.40  0.0
2  80.0  3.05  -0.5
3  36.0  3.31  0.0
```

Entonces `PotencialFF` lee correctamente `NKIND=3, REDELEC=0.35`.
Pero `Potencial.f90` lee `NKIND=3` (ignora el 0.35 del resto de la línea, lo cual en Fortran free-format funciona), así que en este caso **no habría error**.

Sin embargo, si `LJ.dat` tiene el formato que espera `Potencial`:

```
3
1  28.0  3.40  0.0
2  80.0  3.05  -0.5
3  36.0  3.31  0.0
```

Entonces `PotencialFF` intentará leer `NKIND` y `REDELEC` de la primera línea `"3"`, y:
- `NKIND = 3` ✓
- `REDELEC` → falla el read con `IOSTAT /= 0`, y la subrutina retorna sin calcular nada. Las tablas `USS` quedan en cero.

**Consecuencia:** Si `LJ.dat` tiene solo un valor en la primera línea, `PotencialFF` no calcula las tablas de potencial fluido-fluido (`USS`). Todas las energías de interacción entre adsorbatos serán 0, lo que hace que toda molécula insertada sea aceptada (sin repulsión). La simulación pierde todo sentido físico.

Si `LJ.dat` tiene dos valores, la lectura funciona en ambos casos (Fortran ignora el exceso en `Potencial.f90`), pero la variable `REDELEC` de `PotencialFF` no se usa en `Potencial`, lo que es inconsistente.

**Verificación necesaria:** Revisar el archivo `LJ.dat` que se usa actualmente y confirmar el formato. Luego unificar la lectura o separar los formatos.

---

### BUG 3 — ALTO: `Out.f90` no protege contra `N(MOLKIND) = 0`

**Archivo:** `Out.f90`, líneas 48–54

**Qué hace el código:**
Intenta eliminar una molécula aleatoria de tipo `MOLKIND`.

```fortran
MOLKIND = INT(RANF(DUMMY) * NMOLEC) + 1
NTRIAL = N(MOLKIND) - 1

b = RANF(DUMMY)
NLOC = INT(REAL(NTRIAL) * b) + 1
IPULL = LOCATE(NLOC, MOLKIND)
```

**El problema:**
Si `N(MOLKIND) = 0` (no hay moléculas de ese tipo):
- `NTRIAL = -1`
- `NLOC = INT(REAL(-1) * b) + 1` donde b está en [0, 1)
- Si `b = 0.5`: `NLOC = INT(-0.5) + 1 = 0 + 1 = 1` (en algunos compiladores `INT(-0.5) = 0`, en otros `-1`)
- Si `b = 0.99`: `NLOC = INT(-0.99) + 1` → podría ser `0` o incluso negativo

Luego `LOCATE(NLOC, MOLKIND)` con `NLOC = 0` o negativo es un **acceso fuera de rango** del array `LOCATE(5000, 10)`.

**Ejemplo concreto:**

Al inicio de la simulación (ensemble=2), todas las especies tienen `N(I) = 0`. El primer paso MC podría elegir `IJ=2` (out). Se llama a `OUT`, que elige un `MOLKIND` cualquiera, y como `N(MOLKIND) = 0`, se intenta acceder a `LOCATE(0, MOLKIND)` o peor. Dependiendo del compilador, esto puede:
- Dar un segfault
- Leer basura de memoria y continuar con un `IPULL` inválido
- Funcionar "por suerte" si la memoria adyacente tiene un 0

**Corrección necesaria:** Agregar después de `NTRIAL = N(MOLKIND) - 1`:

```fortran
IF (N(MOLKIND) <= 0) RETURN
```

---

### BUG 4 — ALTO: `Main.f90` sobrescribe `Z(NMOLEC)` con actividad del agua

**Archivo:** `Main.f90`, líneas 365–368

```fortran
do INMOLEC = 1, NMOLEC
   Z(INMOLEC) = X(INMOLEC)*P*6.023E-4*((ACEL)**3)*VOL
   Z(NMOLEC)  = 55.55*6.023E-4*((ACEL)**3)*VOL       ! ← ¡dentro del bucle!
end do
```

**El problema:**
En cada iteración del bucle, se sobrescribe `Z(NMOLEC)` con el valor de actividad correspondiente al agua pura (55.55 mol/L). Esto significa:

1. Si `NMOLEC = 2` (por ejemplo, metano + agua):
   - Iteración 1: `Z(1) = X(1)*P*...` (metano) ✓, `Z(2) = 55.55*...` (agua forzada)
   - Iteración 2: `Z(2) = X(2)*P*...` (agua calculada), `Z(2) = 55.55*...` (sobrescrita de nuevo)
   - Resultado: `Z(2)` siempre es 55.55*..., independientemente de `X(2)` y `P`.

2. Si `NMOLEC = 1` (solo una especie):
   - `Z(1) = X(1)*P*...` → inmediatamente se sobrescribe con `Z(1) = 55.55*...`
   - **La actividad de la única especie siempre es la del agua**, sin importar qué especie sea.

**Consecuencia:** En el bloque `NESTADO /= 1`, la actividad de la última especie (o la única) nunca refleja la presión ni la fracción molar correcta. La isoterma estará completamente mal.

**¿Es intencional?** Posiblemente fue pensado para un caso específico donde la última especie es agua y se quiere fijar su actividad. Pero la línea debería estar **fuera del bucle** y solo aplicarse si la última especie es efectivamente agua, con un comentario explicativo.

---

### BUG 5 — MEDIO: `WIJ` no inicializada en `Potin.f90` y `Potout.f90`

**Archivos:** `Potin.f90` línea 101, `Potout.f90` líneas 89 y 121

En ambos archivos:

```fortran
DELTW = DELTW + WIJ
```

`WIJ` se declara como variable local (`REAL :: ... WIJ ...`) pero **nunca se le asigna un valor**. En Fortran, las variables locales no inicializadas contienen basura de memoria.

**¿Por qué no falla actualmente?** Porque `DELTW` no se usa fuera de estas subrutinas (no es un argumento de salida ni una variable de módulo). Pero:
1. Es código peligroso: si alguien activa el virial en el futuro, usará basura.
2. Dependiendo del compilador y flags, una variable no inicializada puede causar señales de NaN que se propagan.

**Contexto:** En el código original, `WIJ` se calculaba a partir de la tabla de potencial. Al refactorizar, parece que se eliminó el cálculo del virial de la tabla pero no la acumulación.

---

### BUG 6 — MEDIO: `Potencial.f90` no inicializa `DELTW` en el bucle interno

**Archivo:** `Potencial.f90`, líneas 68–113

```fortran
DO K = -mat/2, mat/2
   RXI = REAL(K) / REAL(mat)
   DELTV = 0.0       ! ← DELTV sí se inicializa
   VIJ = 0.0
   VIJR = 0.0
   
   DO J = 1, NC
      ...
      DELTW = DELTW + WIJ    ! ← DELTW NO se inicializa a 0 en este punto
   END DO
```

`DELTW` se declara al inicio de la subrutina pero no se pone a cero antes de cada punto de la grilla. Se acumula indefinidamente (o usa basura). No tiene impacto directo porque `DELTW` no se usa después, pero si se descomenta o activa el virial, será un problema.

---

### BUG 7 — BAJO: `Main.f90` llama `print_params()` dos veces

**Archivo:** `Main.f90`, líneas 151–153

```fortran
call read_input('input.txt')    ! ← read_input ya llama a print_params() internamente
call print_params()             ! ← segunda llamada
```

`read_input` (línea 118 de `InputParams.f90`) ya invoca `print_params()`. La segunda llamada en Main imprime el resumen duplicado. También ejecuta dos veces:

```fortran
if (ensemble == 0) dp = 0.0_rk
```

No tiene consecuencia funcional (es idempotente), pero ensucia la salida.

---

### PROBLEMA 1 — ALTO: Generador de números aleatorios de período muy corto

**Archivo:** `Ranf.f90`

```fortran
INTEGER :: L, C, M
PARAMETER (L = 1029, C = 221591, M = 1048576)

SEED = MOD(SEED * L + C, M)
RANF = REAL(SEED) / M
```

**El problema:**
Es un generador lineal congruencial con módulo M = 2^20 = 1,048,576. El período máximo es M = ~10^6 números antes de repetirse.

**Ejemplo concreto:**

Una simulación típica con:
- `isot = 20` puntos de isoterma
- `ijpasos = 1000` pasos de promedio
- `ikpasos = 500` pasos MC
- Cada paso MC llama a `RANF` ~10 veces (selección de movimiento, posición, aceptación, etc.)

Total: 20 × 1000 × 500 × 10 = 100,000,000 llamadas a `RANF`.

Con período ~10^6, la secuencia se repite **~100 veces** durante la simulación. Esto introduce correlaciones sistemáticas: las mismas posiciones de inserción, las mismas decisiones de aceptación/rechazo, etc. Las isotermas pueden mostrar artefactos periódicos.

**Referencia:** El estándar mínimo recomendado para MC es un período de al menos 2^31 (~2×10^9). Idealmente 2^64 o más.

---

### PROBLEMA 2 — MEDIO: Mezcla de precisión simple y doble en todo el código

**Archivos:** Todos

El código trabaja en dos precisiones:
- **Simple (`REAL`):** Energías (V, VG, VA), posiciones (RX, RY, RZ), temperatura, sigma, etc.
- **Doble (`REAL(rk)`):** Módulos PBC (celda, min_image, conversiones frac↔cart).

Las conversiones ocurren en cada paso MC:

```fortran
! En Potin.f90, por cada par de átomos:
dr = min_image(cellR, [real(RX(JIN,K,I),rk), ...], [real(RXI,rk), ...])
RXIJ = real(dr(1))       ! vuelve a simple precisión
```

**Ejemplo de pérdida de precisión:**

```
Valor doble: dr(1) = 0.12345678901234567
Después de real(dr(1)): RXIJ = 0.1234568   (7 dígitos significativos)
```

En el cálculo de energías de Lennard-Jones con `SR6 * (SR6 - 1.0)`, donde `SR6` puede ser ~10^-6, la resta amplifica el error relativo. Para una simulación larga con acumulación de energías, esto puede sesgar los calores isostéricos que dependen de diferencias `<UN> - <U><N>`.

---

### PROBLEMA 3 — MEDIO: Dimensiones hardcodeadas inconsistentes

**Archivos:** Varios

| Constante | Valor | Dónde se usa | Riesgo |
|-----------|-------|-------------|--------|
| `LOCATE(5000, 10)` | 5000 mol, 10 especies | SimulationData | Si NMOLEC > 10, desbordamiento |
| `Z(10)` | 10 especies | In, Out, Move, change | Si NMOLEC > 10, desbordamiento |
| `ANPROM(10)` | 10 especies | Main | Si NMOLEC > 10, desbordamiento |
| `RXBE(50)` | 50 átomos | In, Move, change | Si maxAtoms > 50, desbordamiento |
| `NMAX = 15000` | 15000 moléculas | In.f90 | Inconsistente con 5000 de LOCATE |
| `NMAX = 5000` | 5000 moléculas | Out, Move, change, Add | Consistente con LOCATE |
| `ANX(5000, 10)` | fijo | SimulationData | No es allocatable |

**Ejemplo de fallo:**
Si `NMOLEC = 12` (12 especies de adsorbato), al llamar a `IN`:

```fortran
SUBROUTINE IN(TEMP, Z, SIGMA, EPS, RCUT, ...)
   REAL :: Z(10)     ! ← solo 10 elementos
   ...
   MOLKIND = INT(RANF(DUMMY) * NMOLEC) + 1   ! MOLKIND puede ser 11 o 12
   ...
   DELTCB = BETA * (...) - LOG(Z(MOLKIND) / ...)   ! Z(11) → fuera de rango
```

---

### PROBLEMA 4 — BAJO: `RotationModule` con PI de baja precisión y sin protección de índice

**Archivo:** `RotationModule.f90`

```fortran
REAL, PARAMETER :: PI = 3.14159           ! solo 6 dígitos
REAL, PARAMETER :: PASO = (2.0*PI) / (TABLA_SIZE - 1)
```

El valor correcto de PI a precisión simple es `3.1415927`. Faltan dos dígitos. Además:

```fortran
INDICE_DX = INT((DX + PI) / PASO) + 1
```

Si por errores de redondeo `DX` es ligeramente mayor que `PI` (por ejemplo `3.14160`), entonces:
- `DX + PI = 6.28319`
- `INDICE_DX = INT(6.28319 / 0.006289) + 1 = INT(999.06) + 1 = 1000` ✓

Pero si `DX = PI` exacto (3.14159):
- `DX + PI = 6.28318`
- `INDICE_DX = INT(6.28318 / 0.006289) + 1 = INT(998.9) + 1 = 999` ✓

Y si `DX` es ligeramente mayor que `PI` por redondeo:
- `INDICE_DX` podría ser 1001 → **fuera de rango** de `TABLE_SEN(1000)`.

No hay clamp de seguridad.

---

### PROBLEMA 5 — BAJO: Archivo `Ener###.TXT` se abre y cierra vacío

**Archivo:** `Main.f90`, líneas 389, 452

```fortran
do JPASOS = 1, ijpasos
   ...
   open(unit=51, file='Ener'//CONFIG//'.TXT')    ! abre en cada iteración
   ...
   ! (no se escribe nada en unit 51)
   ...
   close(51)                                      ! cierra vacío
end do
```

El archivo se abre y cierra `ijpasos` veces por isoterma, sin escribir nada. Probablemente quedó de una versión anterior donde se escribían energías parciales.

---

## PARTE 2: PLAN DE CORRECCIÓN

---

### Fase 1 — Bugs críticos (afectan resultados)

Estos deben corregirse antes de ejecutar cualquier simulación de producción.

#### Paso 1.1: Corregir compresión de `EPSAC` y `SGC` en `Estructura.f90`

1. Abrir `Estructura.f90`.
2. Dentro del bloque `IF (all(abs(s) <= 0.5_rk)) THEN` (después de la línea que asigna `QAC(imax) = QAC(i)`), agregar:

```fortran
EPSAC(imax) = EPSAC(i)
SGC(imax)   = SGC(i)
```

3. Verificar que el archivo truncado también se escriba con los valores correctos. Actualmente la línea 114 escribe `EPSAC(i)` y `SGC(i)` con índice `i` del bucle `1..NC` (NC ya actualizado). Esto debería funcionar correctamente después de la corrección.
4. **Test:** Crear un archivo de superficie pequeño (5 átomos) donde 2 estén fuera de la celda. Imprimir `EPSAC(1..NC)` y `SGC(1..NC)` después del filtrado y verificar que correspondan a los átomos dentro de la celda.

#### Paso 1.2: Unificar lectura de `LJ.dat`

1. Determinar el formato correcto de `LJ.dat` revisando archivos de ejemplo existentes.
2. **Opción A** (si `LJ.dat` tiene dos valores en la primera línea): Cambiar `Potencial.f90` línea 35 para que también lea dos valores:

```fortran
READ(11, *) NKIND, REDELEC_dummy    ! o simplemente NKIND (Fortran ignora el resto)
```

En este caso la lectura actual de `Potencial.f90` ya funciona en Fortran free-format (ignora el exceso), pero es mejor ser explícito.

3. **Opción B** (si `LJ.dat` tiene un solo valor en la primera línea): Cambiar `PotencialFF.f90` para no leer `REDELEC` de la primera línea, o leerlo de otra fuente. `REDELEC` se usa como radio de corte electrostático (0.35), pero en `PotencialFF` ya está hardcodeado `RCELE = 0.35`. Se podría eliminar la lectura de `REDELEC`.
4. **Test:** Verificar que ambas subrutinas lean el mismo `NKIND` y los mismos parámetros atómicos.

#### Paso 1.3: Agregar guard en `Out.f90` para `N(MOLKIND) = 0`

1. Abrir `Out.f90`.
2. Después de `MOLKIND = INT(RANF(DUMMY) * NMOLEC) + 1` (línea 48), agregar:

```fortran
IF (N(MOLKIND) <= 0) RETURN
```

3. Verificar que `change.f90` ya tiene protección (sí la tiene: `IF (NTRIAL < 0) RETURN` en línea 82, pero `NLOC` se calcula después, así que está protegido).
4. **Test:** Ejecutar con ensemble=2 (arranque en cero) y verificar que no hay segfaults en los primeros pasos.

#### Paso 1.4: Corregir `Z(NMOLEC)` en `Main.f90`

1. Revisar con el usuario la intención de la línea `Z(NMOLEC) = 55.55*...`.
2. **Si es para agua como última especie:** Mover la línea fuera del bucle y agregar un comentario:

```fortran
do INMOLEC = 1, NMOLEC
   Z(INMOLEC) = X(INMOLEC)*P*6.023E-4*((ACEL)**3)*VOL
end do
! Forzar actividad de la última especie (agua) a concentración pura
Z(NMOLEC) = 55.55*6.023E-4*((ACEL)**3)*VOL
```

3. **Si no es intencional:** Eliminar la línea.
4. **Test:** Imprimir `Z(1:NMOLEC)` y verificar que los valores sean físicamente razonables.

---

### Fase 2 — Bugs de severidad media (código peligroso pero sin impacto inmediato)

#### Paso 2.1: Eliminar acumulación de `WIJ` no inicializado

En `Potin.f90` y `Potout.f90`, eliminar las líneas:

```fortran
DELTW = DELTW + WIJ
```

Ya que `WIJ` no se calcula y `DELTW` no se usa. Si en el futuro se necesita el virial, se debe reimplementar el cálculo completo de `WIJ` a partir de la derivada del potencial.

#### Paso 2.2: Inicializar `DELTW` en `Potencial.f90`

Agregar `DELTW = 0.0` junto con `DELTV = 0.0` dentro del bucle interno, o bien eliminar `DELTW` por completo si el virial no se necesita.

#### Paso 2.3: Eliminar segunda llamada a `print_params()`

En `Main.f90`, eliminar la línea 153:

```fortran
call print_params()    ! ← eliminar, ya se llama desde read_input
```

---

### Fase 3 — Mejoras importantes (calidad de resultados)

#### Paso 3.1: Reemplazar el generador de números aleatorios

Opciones (de menor a mayor esfuerzo):

1. **Mínimo:** Reemplazar el cuerpo de `RANF` por una llamada a `RANDOM_NUMBER` (intrínseco de Fortran 90):

```fortran
REAL FUNCTION RANF(DUMMY)
   IMPLICIT NONE
   REAL :: DUMMY
   CALL RANDOM_NUMBER(RANF)
END FUNCTION RANF
```

Esto usa el generador del compilador (generalmente período ≥ 2^31). Se mantiene la interfaz `RANF(DUMMY)` para no cambiar el resto del código.

2. **Mejor:** Implementar un xoshiro256 o Mersenne Twister en un módulo aparte.
3. **Test:** Ejecutar con el nuevo generador y comparar isotermas. Verificar que el histograma de números generados sea uniforme.

#### Paso 3.2: Hacer dinámicas las dimensiones hardcodeadas

Reemplazar los arrays fijos por allocatables:

- `LOCATE(5000, 10)` → `LOCATE(:,:)` allocatable, dimensionado con `(NMAX_MOL, NMOLEC)`
- `ANX(5000, 10)` → `ANX(:,:)` allocatable
- `ANPROM(10)` → `ANPROM(:)` allocatable con `NMOLEC`
- `Z(10)` en las interfaces → `Z(:)` o `Z(NMOLEC)`
- `RXBE(50)` → `RXBE(:)` allocatable con `maxAtoms`

Esto requiere cambios en las interfaces de `In`, `Out`, `Move`, `change` y en `SimulationData`.

---

### Fase 4 — Mejoras deseables (precisión y robustez)

#### Paso 4.1: Unificar precisión a doble

Cambiar todas las variables de simulación de `REAL` a `REAL(rk)`. Esto es un cambio extenso que toca todos los archivos, pero elimina las conversiones constantes y mejora la precisión de los calores isostéricos.

#### Paso 4.2: Agregar clamp de seguridad en `RotationModule`

```fortran
INDICE_DX = MIN(TABLA_SIZE, MAX(1, INT((DX + PI) / PASO) + 1))
```

#### Paso 4.3: Mejorar PI en `RotationModule`

```fortran
REAL, PARAMETER :: PI = 3.14159265358979
```

O mejor, calcularlo:

```fortran
REAL, PARAMETER :: PI = ACOS(-1.0)
```

#### Paso 4.4: Resolver archivo `Ener###.TXT` vacío

Decidir si se quiere escribir energías por paso (y agregar el write) o eliminar el open/close del unit 51.

---

### Orden de ejecución recomendado

```
Fase 1 (inmediato, antes de cualquier simulación):
  1.1 → 1.2 → 1.3 → 1.4
  Compilar con: bash run2
  Test rápido: ejecutar 1 punto de isoterma con pocos pasos

Fase 2 (limpieza, mismo día):
  2.1 → 2.2 → 2.3
  Compilar y verificar que no hay nuevos warnings

Fase 3 (mejoras, siguiente sesión):
  3.1 (generador aleatorio) → test de regresión
  3.2 (dimensiones dinámicas) → test de regresión

Fase 4 (largo plazo):
  4.1 (precisión doble) → requiere test exhaustivo
  4.2, 4.3, 4.4 → menores, se pueden hacer en cualquier momento
```

---

*Documento generado: febrero 2026. Referencia: lectura completa de todos los archivos .f90 del directorio modulos_f90.*
