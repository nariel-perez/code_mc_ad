# Comparación: logaritmico vs modulos_f90

Este documento compara **explícitamente** los códigos de las carpetas **logaritmico** (Fortran antiguo, archivos `.for`) y **modulos_f90** (Fortran 90/95 modular) del proyecto de simulación Monte Carlo en ensemble gran canónico (GCMC) para adsorción. El objetivo es verificar que ambos realizan, en esencia, la misma simulación.

---

## 1. Flujo y física (equivalente)

| Aspecto | logaritmico | modulos_f90 |
|--------|-------------|-------------|
| **Entrada** | Lee `input.txt` en Main (P, dp, sigmetano, eps, ACEL, T, nam, isot, ijpasos, ikpasos, ensemble, NCELLMAT, NESTADO, etc.) | Igual vía `read_input('input.txt')` |
| **Presiones logarítmicas** | `p_ratio = (dp/p)^(1/(isot-1))`, `p_vals(i) = p_vals(i-1)*p_ratio`, `p = p_vals(ipasos+1)` | Misma fórmula y uso |
| **Unidades reducidas** | SIGMA=sigmetano/ACEL, TEMP=T/eps, PRED=P*SIGMA³/eps, RCUT=10*SIGMA, VOL=XMAX*YMAX*ZMAX | Idéntico |
| **Moléculas** | Lee `MOLEC.DAT`, archivos por molécula, RX0/RY0/RZ0 = (x1/sigmetano)*SIGMA, NATOM, NMIN/NMAXI | Igual vía `read_adsorbates('MOLEC.DAT', ...)` |
| **Superficie y potenciales** | `estructura(eps,nam,sigma,sigmetano,NC,diel)`, `POTENCIALFF(...)`, `POTENCIAL(...)` | Mismas subrutinas y argumentos |
| **Ensemble** | 0=canónico (solo move), 1=GC con initconf, 2=GC desde cero, 3=NAMD | 0, 1, 2 iguales; 3 eliminado en modulos_f90 |
| **Restart** | Si ensemble<2: lee `initconf.txt` (V,VG,VA, nmolec2, por tipo natom2/ncantmol, coords, `call add(molkind)`) | Misma lógica |
| **Energía inicial** | Bucle sobre moléculas: `adpotout` + `potout` | Mismo bucle y llamadas |

---

## 2. Bucle principal (isotermas → promedios → pasos MC)

- **Isotermas**: `IPASOS=1,isot`; reset de CNF; apertura de `CONFIG###.xyz` y `CONFIG###.TXT`; cálculo de actividades **Z** (NESTADO 1 vs otro); puesta a cero de acumuladores (U, UG, UA, UN, etc.).
- **Promedios**: `JPASOS=1,ijpasos`; reescalado V/VG/VA a unidades reducidas; `MULT=mult2` solo si JPASOS==1.
- **Pasos MC**: `KPASOS=1, ikpasos*MULT`.
  - Si **ensemble==0**: solo `call move(...)`.
  - Si no: elección aleatoria del movimiento:
    - **logaritmico**: `ij = ranf(dummy)*3+1` → solo 1, 2 o 3; luego `goto (10,20,30,35) ij` (10=in, 20=out, 30=move, 35=change). El caso 4 (change) casi no se elige porque el aleatorio es *3+1.
    - **modulos_f90**: `IJ = int(ranf(DUMMY)*4)+1` y `select case (1:in, 2:out, 3:move, 4:change)` → cuatro movimientos con probabilidad correcta.
- Tras KPASOS: reescalado V/VG/VA; acumulación de `NTOTAL(I)`, `estadistica(CNF)`; acumuladores U, UG, UA, UN, UNG, UNA, u2, AN, N2. Misma lógica en ambos.

---

## 3. Subrutinas de movimiento (misma función)

- **In**: elegir tipo de molécula (MOLKIND), posición aleatoria (RXNEW,RYNEW,RZNEW), cargar coordenadas base (RX0→RXBE/RX1), rotación (en logaritmico con ángulos/COS/SIN, en modulos_f90 con `GetRotationMatrix`), **potin** (DELTV), **adpotin** (DELTVA), criterio Metropolis, si se acepta **add(molkind)**. Misma secuencia.
- **Out**: elegir molécula, **potout** + **adpotout**, Metropolis, si se acepta **remove**. Igual.
- **Move** y **change**: misma idea (potout+potin, adpotout+adpotin, aceptación, actualización de posiciones o tipo).

Las interfaces pueden diferir (p. ej. en logaritmico `Out` recibe `NMIN,NMAXI`; en modulos_f90 se pasan menos argumentos por módulos). La lógica de “calcular cambio de energía y aceptar/rechazar” es la misma.

---

## 4. Salidas por isoterma

- Escritura de **initconf.txt** (V,VG,VA y configuración por tipo).
- **Promedios**: `ANPROM(I) = NTOTAL(I)/ijpasos` (en logaritmico se usa `JPASOS` en el denominador, que es el contador del bucle, equivalente a `ijpasos`).
- **Calores isostéricos**: mismas fórmulas (AN1, U1, UNG1, UG1, UNA1, UA1, UN1, u2, AN2, ANN, CALOR, CALORG, CALORA, CALORESP1/2/3).
- **XYZ**: lectura del archivo de superficie truncada (átomos de la estructura + moléculas adsorbidas). En logaritmico el archivo se llama **`TRUNCADO.TXT`** (nombre fijo); en modulos_f90 **`<base>_truncado.txt`** (derivado de `nam`). En ambos, Estructura escribe ese archivo y Main lo lee para armar el XYZ.
- **CNF**: escritura de `CNF###-I-NATOMKINDI.TXT` con la misma grilla (NCELLMAT).
- **SALIDAACTIVADO-100.TXT**: P, ANPROM, CALOR, CALORA, CALORESP3.
- Actualización de presión: `p = p_vals(ipasos+1)` en ambos.

---

## 5. Resumen de diferencias (sin cambiar la esencia)

| Aspecto | logaritmico | modulos_f90 |
|--------|-------------|-------------|
| **Organización** | COMMON blocks (BLOCK1, BLOCK2, BLOCK3) | Módulos (InputParams, SimulationData, AdsorbateInput, etc.) |
| **Lectura** | Todo en Main | `read_input` y `read_adsorbates` |
| **Elección de movimiento** | `*3+1` → “change” casi no se usa | `*4+1` → cuatro movimientos equiprobables |
| **Control de flujo** | `goto (10,20,30,35)` | `select case` |
| **Archivo truncado** | Nombre fijo `TRUNCADO.TXT` | `<base>_truncado.txt` |
| **Código NAMD** | ensemble=3 y NAMD1/NAMD2 | Eliminado |
| **Cierres de archivo** | `close(21)`, `close(22)` (unidades no abiertas en el flujo mostrado) | Limpiado |

---

## 6. Conclusión

Los códigos de **logaritmico** y **modulos_f90** implementan, en esencia, la **misma simulación GCMC de adsorción**: mismas entradas, mismas ecuaciones (unidades reducidas, actividades, calores isostéricos), mismos movimientos MC (in, out, move, change) y mismas salidas (isotermas, configuraciones, CNF, initconf). La versión en **modulos_f90** es una modernización (módulos, lectura en subrutinas, select case, nombre dinámico del truncado y corrección del muestreo de los cuatro movimientos), sin cambiar la física ni el algoritmo central.

---

## 7. Diferencias numéricas en la salida (SALIDAACTIVADO-100.TXT)

Al correr el **mismo sistema** (p. ej. adsorción de butano en grafito) con ambos códigos, las filas de salida pueden diferir. Esta sección explica a qué se debe y por qué las “energías” (columnas 3–5) pueden cambiar mucho aunque la diferencia en número de moléculas (columna 2) sea pequeña.

### 7.1 Formato de la salida

Cada fila de `SALIDAACTIVADO-100.TXT` corresponde a un punto de la isoterma y tiene la forma:

`P` | `ANPROM(1)` | `CALOR` | `CALORA` | `CALORESP3`

- **Columna 1**: presión (P).
- **Columna 2**: número promedio de moléculas adsorbidas, ANPROM(1).
- **Columnas 3–5**: magnitudes derivadas de promedios y fluctuaciones de energía y N:
  - **CALOR**: calor isostérico (total).
  - **CALORA**: contribución adsorbato–superficie al calor isostérico.
  - **CALORESP3**: cantidad relacionada con fluctuaciones (CALORESP1/CALORESP2).

### 7.2 Por qué puede diferir el número de moléculas (columna 2)

- **Rotación**  
  - **logaritmico**: usa ángulos (DX, DY, DZ) y calcula la matriz de rotación con **COS/SIN directos** (continuos).  
  - **modulos_f90**: usa **tabla precalculada** de 1000 valores de seno/coseno en `[−π, π]`; el ángulo se discretiza al más cercano (`INDICE = INT((DX+PI)/PASO)+1`).  
  Eso introduce un pequeño error en la matriz de rotación → posiciones de prueba ligeramente distintas → aceptaciones/rechazos distintos → trayectoria MC distinta y, en consecuencia, un **N promedio** que puede diferir en unas pocas moléculas (p. ej. ~3 al final de la isoterma).

- **Elección del movimiento**  
  En logaritmico el movimiento se elige con `*3+1` (solo 1, 2 o 3), por lo que “change” casi no se usa; en modulos_f90 se usa `*4+1` (in, out, move, change equiprobables). Eso también hace que las trayectorias no sean las mismas.

Ambos efectos son de **implementación** (rotación tabulada vs continua, y mezcla de movimientos), no de la física del modelo. Que ANPROM coincida bastante (sobre todo a presiones bajas) y solo se separe un poco a alta presión es coherente con eso.

### 7.3 Por qué CALOR, CALORA y CALORESP3 pueden cambiar mucho

Esas tres columnas **no son** la energía total en sí, sino cantidades derivadas de **promedios y fluctuaciones** a lo largo del run:

- **CALOR** = 8.3144*T − ((UN1 − U1*AN1) / **ANN**)*8.31  
- **CALORA** = −((UNA1 − UA1*AN1) / **ANN**)*8.31  
- **CALORESP1** = (u2 − U1²)*8.31² − (8.31*U1*AN1)² / **ANN**  
- **CALORESP2** = UN1 − AN1*8.31*T²  
- **CALORESP3** = CALORESP1 / CALORESP2  

Donde **ANN = AN2 − AN1²** es la **varianza del número de partículas** (fluctuación de N). Aparece en el **denominador** de CALOR y CALORA. Por tanto:

1. **Trayectorias distintas**  
   Aunque se use la misma semilla de números aleatorios en ambos códigos, la rotación tabulada y la distinta mezcla de movimientos hacen que, desde el primer paso MC, las configuraciones aceptadas sean distintas. Es decir, son **dos realizaciones distintas** del mismo proceso estocástico. Entonces U, UA, UG, UN, UNA, etc., y sobre todo sus **fluctuaciones** (u2, AN2, covarianzas), no tienen por qué ser similares entre sí.

2. **Sensibilidad a las fluctuaciones**  
   CALOR y CALORA dependen de cocientes donde entra **ANN** (varianza de N). Si ANN es pequeña, pequeños cambios en UN1, U1, AN1, UNA1, UA1 (por tener otra trayectoria) producen **grandes cambios** en (UN1−U1*AN1)/ANN y (UNA1−UA1*AN1)/ANN. Lo mismo para CALORESP1 (que también divide por ANN) y para CALORESP3. Por eso es normal que, **incluso con una diferencia pequeña en ANPROM** (p. ej. 3 moléculas), las columnas 3–5 difieran bastante: no solo cambia un poco el promedio de N, sino toda la estadística de fluctuaciones (U, N, productos U×N, etc.).

3. **Resumen**  
   Las “energías” que ves en las columnas 3–5 son **calores isostéricos y cantidades relacionadas con fluctuaciones**. Dependen mucho de la trayectoria MC concreta y de varianzas/covarianzas. Que cambien mucho entre logaritmico y modulos_f90, aunque ANPROM sea parecido, es **esperable** por:  
   - trayectorias diferentes (rotación + mezcla de movimientos),  
   - fórmulas con ANN en el denominador (muy sensibles a pequeñas diferencias en promedios y fluctuaciones).

Si se quisiera acercar más los resultados entre ambos códigos, habría que: (i) usar la misma rotación (seno/coseno continuos también en modulos_f90, sin tabla), y (ii) igualar la elección de movimientos (*3+1 y mismo orden de casos que en logaritmico, o aceptar que siempre habrá diferencias por ser dos realizaciones estocásticas distintas).

---

*Documento generado a partir de la revisión comparativa de los directorios `logaritmico/` y `modulos_f90/` del proyecto code_mc_ad.*
