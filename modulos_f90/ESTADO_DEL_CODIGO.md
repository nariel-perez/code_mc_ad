# Estado del código — modulos_f90

Documento de referencia con el **trabajo realizado** y el **estado actual** del proyecto de simulación Monte Carlo para adsorción (GCMC).

---

## 1. Descripción del proyecto

- **Objetivo**: Simulación de Monte Carlo en ensemble gran canónico (GCMC) para adsorción de gases en superficies porosas (p. ej. carbón activado).
- **Lenguaje**: Fortran 90/95 (módulos, subrutinas).
- **Entradas**: `input.txt`, `MOLEC.DAT`, archivo de superficie (nombre en `nam`), opcionalmente `initconf.txt` para restart.
- **Salidas**: Isotermas, calores isostéricos, configuraciones XYZ, estadísticas espaciales (CNF), archivo truncado de superficie, `initconf.txt` para restart.

---

## 2. Trabajo realizado (resumen)

### 2.1 Mejoras en Main.f90 (registradas en MEJORAS.md)

- Corrección de llamada a `out()` (eliminación de parámetros `nmin`, `nmaxi`).
- Índice `I` en bucle para `anprom`/`ntotal`.
- Eliminación de `close(21)` y `close(22)` (archivos no abiertos).
- Eliminación de variables no usadas, redundancias (`P = P`, `MULT`, `escalax`/`escalay`/`escalaz`, `ESCALA`, bloque vacío `ensemble2`).
- Sustitución de `goto` por `select case` y corrección del bug en la elección de movimiento (`*3` → `*4` para cuatro casos: in, out, move, change).
- Eliminación de `write(58,*)` sin unidad abierta y de variables `aitest76`/`aitest77` (comparación con real sustituida por `mod(JPASOS, 500) == 0`).
- Inicialización de `u2` y conversión explícita de `p_ratio`; eliminación de truncado de caracteres en CONFIG/CONFAT/CONFNAT.
- Eliminación de variables locales no usadas (p. ej. `RXNEW`, `RYNEW`, `RZNEW`, `DELTW` en Main).

### 2.2 Bugs y conversiones (varios archivos)

- **In.f90**: Inicialización de `DELTW = 0.0` para evitar uso no inicializado.
- **Main.f90**: Inicialización de `u2` en el bloque de acumuladores por isoterma.
- **Potin.f90, Potout.f90**: Conversión explícita REAL→INTEGER para índice de tabla: `IDIST = INT(RIJ*1000.0 + 1.0)`.

### 2.3 Archivo truncado (Estructura + Main)

- **Estructura.f90**: El archivo de superficie truncada se escribe con nombre `<base>_truncado.txt`, donde `base` es el nombre del archivo de superficie sin extensión (ej. `carbo.txt` → `carbo_truncado.txt`). Se evita doble extensión.
- **Main.f90**: Variable `archivo_truncado` construida con la misma regla; apertura con `open(unit=49, file=archivo_truncado)` para leer la superficie truncada al escribir configuraciones XYZ.

### 2.4 Código obsoleto eliminado

- Eliminado el bloque `if (ensemble.eq.3) then call namd1 ... ensemble = 1 end if` y la variable `ensemble2` en Main (las rutinas `namd1`/`namd2` ya no existen).
- Eliminadas las líneas comentadas de compilación de `namd1.f` y `namd2.f` en el script `run2`.

### 2.5 Limpieza de warnings

- Eliminadas variables y parámetros no usados en múltiples módulos (PotencialFF, Potencial, In, Adpotin, Potin, Add, Out, Potout, Adpotout, Remove, Move, change), manteniendo las variables necesarias en el cuerpo del código.
- No se usan directivas `!GCC$ ATTRIBUTES UNUSED` porque el compilador actual no las reconoce; persisten avisos de argumentos dummy no usados (sin impacto en ejecución).

### 2.6 Documentación y compilación

- Unificación de la documentación de mejoras en un único archivo: `MEJORAS.md` (eliminados `MEJORAS_SUGERIDAS.md` y `ANALISIS_MEJORAS_PENDIENTES.md`).
- **run2**: Al final de la compilación se ejecuta `rm -f *.o *.mod` para no dejar restos de compilación (solo fuentes y, si se desea, el ejecutable).

---

## 3. Estado actual del código

### 3.1 Compilación

- **Comando**: `bash run2` (desde el directorio del proyecto).
- **Resultado**: Compila sin errores; genera `simulation.exe`.
- **Avisos**: Warnings por argumentos dummy y algunas variables no usadas (principalmente en subrutinas In, Out, Move, change y en módulos de potencial). No afectan la corrección de los resultados.

### 3.2 Estructura de archivos

| Tipo        | Archivos |
|------------|----------|
| **Módulos** | InputParams.f90, AdsorbateInput.f90, SimulationData.f90, ConstantsModule.f90 (PhysicalConstants), RotationModule.f90, Estructura.f90 (EstructuraModule) |
| **Potenciales** | PotencialFF.f90, Potencial.f90 |
| **Movimiento MC** | In.f90, Out.f90, Move.f90, change.f90 |
| **Auxiliares** | Add.f90, Remove.f90, Adpotin.f90, Adpotout.f90, Potin.f90, Potout.f90, estadistica.f90, Ranf.f90 |
| **Principal** | Main.f90 |
| **Build** | run2 |
| **Documentación** | DIAGRAMA_FLUJO.md, flujo_claude.md, MEJORAS.md, ESTADO_DEL_CODIGO.md; opcional: compilation_warnings.txt |

### 3.3 Flujo resumido

1. Lectura de entrada (`input.txt`, `MOLEC.DAT`, archivo de superficie).
2. Constantes físicas y tablas de rotación.
3. Lectura de adsorbatos y alocación de arrays.
4. Estructura de superficie → escribe `<base>_truncado.txt`; Main guarda el nombre en `archivo_truncado`.
5. Construcción de tablas de potencial (PotencialFF, Potencial).
6. Inicialización según ensemble (arranque en cero o desde `initconf.txt`).
7. Bucle de isotermas (IPASOS) → bucle de promedios (JPASOS) → bucle de pasos MC (KPASOS). En cada paso: si `ensemble == 0` solo move; si no, `select case (IJ)` con 1=in, 2=out, 3=move, 4=change.
8. Estadísticas, calores isostéricos, escritura de configuraciones y CNF; actualización de presión; guardado de `initconf.txt`.
9. Cierre de archivos y desalocación.

### 3.4 Archivos de entrada/salida (actual)

- **Entrada**: `input.txt`, `MOLEC.DAT`, archivos de moléculas, archivo de superficie (`nam`), `initconf.txt` (si ensemble &lt; 2).
- **Salida**: `<base>_truncado.txt`, SALIDAACTIVADO-100.TXT, PERFILES.TXT, CONFIG###.xyz, CONFIG###.TXT, CNF###-*-*.TXT, Ener###.TXT, initconf.txt.

---

## 4. Cómo compilar y ejecutar

```bash
cd /ruta/modulos_f90
bash run2
./simulation.exe   # con input.txt, MOLEC.DAT y archivo de superficie en su lugar
```

Tras compilar, `run2` borra `.o` y `.mod`; solo quedan los fuentes y `simulation.exe` (o solo fuentes si se elimina el ejecutable).

---

## 5. Referencia a la documentación

| Documento | Contenido |
|-----------|-----------|
| **DIAGRAMA_FLUJO.md** | Diagrama detallado del flujo, módulos, subrutinas MC y archivos de E/S. |
| **flujo_claude.md** | Resumen del flujo GCMC, bloques principales y diagrama simplificado. |
| **MEJORAS.md** | Registro de mejoras aplicadas (Main y resto), bugs corregidos y trabajo posterior. |
| **ESTADO_DEL_CODIGO.md** | Este archivo: trabajo realizado y estado actual del código. |

---

*Última actualización: febrero 2025.*
