# MEJORAS EN Main.f90 — Registro único

Documento unificado de mejoras aplicadas al programa principal. Sirve como referencia para revisiones futuras.

---

## Estado actual

**Todas las mejoras listadas han sido aplicadas** en `Main.f90`. El código compila y no presenta los errores ni redundancias descritos a continuación (que existían en versiones anteriores).

---

## 1. Problemas críticos (corregidos)

### 1.1 Parámetros incorrectos en `call out()`
- **Problema**: Se pasaban `nmin` y `nmaxi` a `OUT`, pero la subrutina no los acepta (firma en `Out.f90`).
- **Solución aplicada**: Llamada con `(temp, z, sigma, eps, rcut, v, va, vg, w, ghost, jpasos, canonicalmolecules)`.

### 1.2 Índice fuera de bucle (anprom / ntotal)
- **Problema**: `anprom(I)` y `ntotal(I)` se asignaban fuera de un bucle, usando un `I` indefinido.
- **Solución aplicada**: Encerrado en `do I = 1, NMOLEC` con `anprom(I) = 0` y `ntotal(I) = 0`.

### 1.3 Cierre de archivos no abiertos
- **Problema**: `close(21)` y `close(22)` sin correspondiente `open` en Main.
- **Solución aplicada**: Eliminadas esas líneas.

---

## 2. Variables no utilizadas (eliminadas)

Se eliminaron las declaraciones locales no usadas en Main:  
`NMAX`, `RMIN`, `OVRLAP`, `MOLEC1`, `NMATOM`, `IKIND`, `NS`, `NCONFMIN`, `NCONFMAX`.

**Nota**: En una revisión posterior se eliminaron también `RXNEW`, `RYNEW`, `RZNEW` y `DELTW` de Main (no se usan en el cuerpo del programa).

---

## 3. Redundancias y limpieza (aplicadas)

| Mejora | Acción realizada |
|--------|-------------------|
| Asignación `P = P ! * 1333.22` | Línea eliminada. |
| Doble asignación de `MULT = 1` | Dejada una sola asignación antes del `if (JPASOS.eq.1)`. |
| `ntotalGRAF = NC` dentro del bucle | Eliminada la asignación redundante; `ntotalGRAF` se usa correctamente (NC + acumulación por moléculas). |
| Variables `escalax`, `escalay`, `escalaz` | Sustituidas por `write(*,*) acelx, acely, acelz`. |
| Variable `ESCALA = acel` | Sustituida por uso directo de `acel` en las expresiones (p. ej. `RX(...)*acel`). |
| Bloque vacío `if (ensemble2.eq.3)` | Bloque eliminado. |

---

## 4. Estilo y organización (aplicadas)

| Mejora | Acción realizada |
|--------|-------------------|
| Líneas en blanco excesivas | Reducidas entre declaraciones e inicio del programa. |
| Bucle `do  i = 2, isot` | Modernizado a `do i = 2, isot`. |
| `goto` con etiquetas numéricas | Sustituido por `if (ensemble == 0)` + `select case (IJ)` con casos 1–4. |
| Bug en elección de movimiento | Corregido: antes `ranf(DUMMY)*3+1` (solo 1–3) con 4 casos; ahora `int(ranf(DUMMY)*4)+1` (1–4). |
| `write(58,*)` sin unidad abierta | Línea eliminada (evita error en ejecución). |

---

## 5. Decisiones mantenidas (sin cambio)

- **Archivo unidad 51** (`Ener...TXT`): Sigue abriéndose y cerrándose dentro del bucle JPASOS; se asume intencional (p. ej. sobrescritura por iteración).
- **Desalocación**: Se mantiene la desalocación actual de arrays (Z, X, N, NATOM, etc.) según el diseño del programa.

---

## 6. Resumen numérico

- Mejoras aplicadas: 14.
- Líneas modificadas o eliminadas: ~25–30.
- Bugs corregidos: 2 (elección de movimiento 1–4; escritura en unidad 58 no abierta).
- Errores de compilación conocidos: 0.
- Referencia de código: `Main.f90` actual.

---

## 7. Trabajo posterior (revisión 2025)

### Bugs y conversiones (corregidos)
- **In.f90**: Inicialización de `DELTW = 0.0` para evitar uso no inicializado en `W = W + DELTW`.
- **Main.f90**: Inicialización de `u2 = 0` en el bloque de acumuladores por isoterma; eliminación de comparación con real (`aitest77.eq.0` → `mod(JPASOS, 500) == 0`); conversión explícita de `p_ratio`; eliminación de truncado de caracteres (CONFIG/CONFAT/CONFNAT); eliminación de variables `aitest76`, `aitest77`.
- **Potin.f90, Potout.f90**: Conversión explícita REAL→INTEGER para índice de tabla: `IDIST = INT(RIJ*1000.0 + 1.0)`.

### Archivo truncado (estructura + Main)
- **Estructura.f90**: El archivo de superficie truncada se escribe con nombre `<base>_truncado.txt`, donde base es `nam` sin extensión (ej. `carbo.txt` → `carbo_truncado.txt`). Se evita doble extensión.
- **Main.f90**: Variable `archivo_truncado` y misma regla de nombre; apertura con `open(unit=49, file=archivo_truncado)`.

### Código obsoleto eliminado
- Eliminado bloque `if (ensemble.eq.3) then call namd1 ... end if` y variable `ensemble2` en Main; eliminadas referencias a `namd1`/`namd2` en el script de compilación.

### Limpieza de warnings
- Eliminadas variables y parámetros no usados en varios módulos (PotencialFF, Potencial, In, Adpotin, Potin, Add, Out, Potout, Adpotout, Remove, Move, change). Se mantuvieron las variables necesarias para el cuerpo del código. Las directivas `!GCC$ ATTRIBUTES UNUSED` no son reconocidas por el compilador actual, por lo que persisten avisos de argumentos dummy no usados.

### Script de compilación (run2)
- Al final de la compilación se ejecuta `rm -f *.o *.mod` para dejar solo fuentes y ejecutable (o solo fuentes si se borra el .exe).

---

*Última actualización: febrero 2025. Para flujo del programa y módulos, ver `DIAGRAMA_FLUJO.md` y `flujo_claude.md`. Para estado actual detallado, ver `ESTADO_DEL_CODIGO.md`.*
