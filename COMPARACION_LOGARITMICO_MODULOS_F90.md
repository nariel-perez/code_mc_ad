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

*Documento generado a partir de la revisión comparativa de los directorios `logaritmico/` y `modulos_f90/` del proyecto code_mc_ad.*
