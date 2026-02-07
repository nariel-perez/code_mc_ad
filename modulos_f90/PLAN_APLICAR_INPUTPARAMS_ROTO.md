# Plan: aplicar estructura de roto/InputParams en modulos_f90/InputParams

Este documento describe **cómo** adaptar el `InputParams.f90` de modulos_f90 para que tenga la misma estructura y capacidades que el de roto (input clave–valor, celda PBC, ángulos, dim), **sin aplicar aún los cambios**. Sirve como guía para la Fase 1 (y preparación Fase 2) del plan de trabajo PBC/input.

---

## 1. Diferencias actuales (resumen)

| Aspecto | modulos_f90 (actual) | roto |
|--------|----------------------|------|
| **Lectura** | Orden fijo: una línea por parámetro en orden (P, dp, sigmetano, …) | Clave–valor: líneas `clave valor` o `clave valor1 valor2`, comentarios `!`, líneas vacías ignoradas |
| **Precisión** | `real` (simple) | `real(rk)` con `rk = selected_real_kind(12,99)` vía `use PBC_Mod` |
| **PBC / celda** | Solo BCX, BCY, BCZ (enteros); no hay tipo Cell | BCX, BCY, BCZ + `type(Cell_t) :: cell`; se construye con `cell_from_lengths_angles` tras leer |
| **Geometría** | acelx, acely, acelz; ACEL | Igual + angux/anguy/anguz (α,β,γ en grados), opcionales; dim_flag (2/3); pointers a_len=>acelx etc. |
| **Rutina de impresión** | `print_params()` | `print_summary()` (incluye ángulos, dim, mensajes estilo roto) |
| **Dependencias** | Ninguna de otros módulos de proyecto | `use PBC_Mod, only : Cell_t, cell_from_lengths_angles, rk` |
| **Auxiliares** | No hay | `lowercase()`, `set_int()`, `set_real()`, `set_char()`, `info()`, `warn()`, `fatal()` |

---

## 2. Pasos concretos para aplicar (orden sugerido)

### 2.1 Añadir dependencia de PBC_Mod

- Al inicio del módulo, después de `module InputParams`:
  - Añadir:  
    `use PBC_Mod, only : Cell_t => Cell, cell_from_lengths_angles, rk`
  - Mantener `implicit none` y, si se desea, `private` con `public` explícito solo para lo que el resto del código usa.

- **Importante**: El orden de compilación debe ser PBC_Mod antes que InputParams (ya previsto en el plan: run2 debe compilar PBC_Mod.f90 antes de InputParams.f90).

### 2.2 Cambiar tipos de los parámetros públicos numéricos

- Sustituir `real` por `real(rk)` en todas las variables públicas reales:  
  P, dp, sigmetano, eps, ACEL, acelx, acely, acelz, diel, T.
- Mantener los nombres que usa el resto del código (Main, In, Move, Potin, Potout, etc.):  
  `acelx`, `acely`, `acelz`, `ACEL`, y en el resto de módulos se usa `bcx`/`bcy`/`bcz` en el `USE` pero en roto son `BCX`/`BCY`/`BCZ`; en Fortran es lo mismo (case-insensitive). Se puede mantener `BCX`, `BCY`, `BCZ` para coincidir con roto y con el uso en el código (In.f90, etc., usan BCX, BCY, BCZ).

### 2.3 Añadir variables de geometría y celda (como en roto)

- Añadir reales para ángulos (por defecto 90°):  
  `real(rk), public, target :: angux=90._rk, anguy=90._rk, anguz=90._rk`
- Opcional (roto lo tiene): pointers para compatibilidad con nombres de PBC_Mod:  
  `real(rk), public, pointer :: a_len=>acelx, b_len=>acely, c_len=>acelz`  
  `real(rk), public, pointer :: alpha_deg=>angux, beta_deg=>anguy, gamma_deg=>anguz`
- Añadir dimensionalidad:  
  `integer, public :: dim_flag = 3`
- Añadir la celda global:  
  `type(Cell_t), public :: cell`

### 2.4 Sustituir el cuerpo de read_input

- Eliminar la lectura por orden fijo (los `read(iu,*) P`, `read(iu,*) dp`, …, `read(iu,*) BCX, BCY, BCZ`, etc.).
- Implementar lectura por líneas:
  1. Abrir con `open(newunit=iu, file=fname, status='old', action='read')` (y manejo de error si falla).
  2. Bucle: `read(iu,'(A)',iostat=ios) line` hasta `iostat_end`.
  3. Para cada línea: quitar comentarios (todo desde `!`), ignorar si `trim(line)==''`.
  4. Intentar leer `key` y un valor: primero como entero (`read(line,*,iostat=ios) key, ival`); si ios==0 → `set_int(lowercase(trim(key)), ival)`; si no, intentar como real → `set_real(lowercase(trim(key)), rval)`; si no, intentar como cadena → `set_char(line)`.
  5. Cerrar el archivo.

- Usar `iso_fortran_env` para `iostat_end` (y opcionalmente `output_unit` en print_summary).

### 2.5 Después de leer todas las claves: ACEL y construcción de cell

- Si `ACEL == 0._rk`, asignar:  
  `ACEL = max(acelx, acely, acelz)` (o max(a_len, b_len, c_len) si se usan los pointers).
- Llamar:  
  `call cell_from_lengths_angles(cell, acelx, acely, acelz, angux, anguy, anguz, dim=dim_flag, centered=.true.)`  
  (usar los nombres de variables que se definan: a_len/b_len/c_len si existen como alias, si no acelx/acely/acelz y alpha_deg/beta_deg/gamma_deg o angux/anguy/anguz).
- Aplicar PBC y dimensión a la celda:  
  `cell%pbc = [ BCX/=0, BCY/=0, BCZ/=0 ]`  
  `cell%dim = dim_flag`
- Llamar a la rutina de resumen:  
  `call print_summary()` (ver 2.7).

### 2.6 Implementar set_int, set_real, set_char y lowercase

- **set_int(key, ival)**: según `key` en minúsculas asignar a: mat, isot, ijpasos, ikpasos, mult2, ensemble, NESTADO, canonicalmolecules, NCELLMAT, BCX, BCY, BCZ, dim_flag. Claves sugeridas (como en roto): 'mat', 'isot', 'ijpasos', 'ikpasos', 'mult2', 'ensemble', 'nestado', 'canonicalmolecules', 'ncellmat', 'bcx', 'bcy', 'bcz', 'dim'. Cualquier otra clave → `call warn('Clave entera desconocida: '//trim(k))`.
- **set_real(key, rval)**: según key asignar a P, dp, sigmetano, eps, diel, T, acelx, acely, acelz, angux, anguy, anguz, ACEL. Claves sugeridas: 'p', 'dp', 'sigmetano', 'eps', 'diel', 't', 'acelx'/'a'/'ax', 'acely'/'b'/'ay', 'acelz'/'c'/'az', 'angux'/'alpha'/…, 'acel'. Cualquier otra → warn.
- **set_char(line)**: leer de la línea key y valor (cadena); si key es 'nam' o 'nombre' asignar a `nam`; si no, warn.
- **lowercase(str)**: función pura que devuelve la cadena con letras mayúsculas convertidas a minúsculas (iachar/achar, rango 'A'–'Z' → +32).

Todas estas subrutinas/funciones son internas al módulo (no hace falta que sean públicas si solo las usa read_input).

### 2.7 Renombrar y ampliar print_params → print_summary

- Renombrar `print_params` a `print_summary` (y que sea la única rutina de resumen que se llame al final de read_input).
- Añadir en la impresión:
  - Ángulos (angux, anguy, anguz o alpha_deg, beta_deg, gamma_deg).
  - dim_flag (dimensionalidad).
- Ajustar textos para que coincidan con el estilo de roto si se desea (por ejemplo "Cell size (Ref,x,y,z):", "Angles (α,β,γ) [deg]:", "BC (Ref,BCx,BCy,BCz):", "dimensionalidad:").
- Mantener la validación actual de NCELLMAT > 50 (stop con mensaje) y la línea que hace `if (ensemble == 0) dp = 0.0`.
- Usar `output_unit` de `iso_fortran_env` para la unidad de escritura si se quiere igualar a roto.

### 2.8 Auxiliares info, warn, fatal

- Añadir tres subrutinas internas (privadas):
  - `subroutine info(msg)`  → `write(*,'(A)') trim(msg)`
  - `subroutine warn(msg)`   → `write(*,'("WARNING: ",A)') trim(msg)`
  - `subroutine fatal(msg)`  → `write(*,'("FATAL: ",A)') trim(msg); stop`
- En read_input: al inicio se puede llamar `call info('Leyendo '//trim(fname))`; si falla la lectura de una línea, `call fatal('Error leyendo línea')` (o equivalente); en set_* se usa `warn` para claves desconocidas.

### 2.9 Compatibilidad con input por orden fijo (opcional, Fase 1)

- Según el plan hay dos opciones:  
  **A)** Detección: si la primera línea leída parece un número (o la primera palabra es numérica), usar lectura por orden fijo; si no, usar clave–valor.  
  **B)** Solo clave–valor y documentar el nuevo formato (con script o nota para convertir inputs viejos).
- Para este plan de “cómo aplicar”: se deja decidido que se implementará **solo clave–valor** como en roto; si más adelante se quiere compatibilidad legacy, se añadirá un `if` al inicio del bucle de read_input que detecte formato y bifurque a una subrutina `read_input_legacy()` que reproduzca los read por orden actuales.

### 2.10 Nombres públicos que debe seguir exportando InputParams

Para no romper el resto de modulos_f90, el módulo debe seguir exportando (y con el mismo nombre):

- Reales: P, dp, sigmetano, eps, ACEL, acelx, acely, acelz, diel, T  
- Enteros: BCX, BCY, BCZ, mat, isot, ijpasos, ikpasos, mult2, ensemble, NESTADO, canonicalmolecules, NCELLMAT  
- Cadena: nam  
- Nuevos: angux, anguy, anguz (o alpha_deg, beta_deg, gamma_deg), dim_flag, cell  
- Subrutinas públicas: read_input, print_summary (y si Main llama a print_params, cambiar en Main la llamada a print_summary o mantener un alias print_params => print_summary).

### 2.11 Cambios en Main.f90 (cuando se aplique)

- Donde se llame `print_params()`, cambiar a `print_summary()` (o añadir en InputParams `print_params` como alias público que llame a print_summary).
- Asegurar que Main (y cualquier otro que use InputParams) no asume `real` de precisión simple para las variables que pasen a `real(rk)` (en la práctica, si todo usa el mismo módulo y rk es el real de doble precisión típico, no suele haber problema).

### 2.12 Orden de compilación en run2

- Incluir en run2, en este orden:  
  1) PBC_Mod.f90  
  2) GeomUtils.f90  
  3) InputParams.f90  
  4) … resto como hasta ahora (ConstantsModule, AdsorbateInput, etc.).

---

## 3. Resumen de archivos a tocar (cuando se haga)

| Archivo | Cambios |
|---------|--------|
| **modulos_f90/InputParams.f90** | use PBC_Mod; real→real(rk); añadir angux/anguy/anguz, dim_flag, cell; reemplazar lectura por clave–valor; set_int/set_real/set_char, lowercase; ACEL por defecto; cell_from_lengths_angles y cell%pbc/dim; print_summary; info/warn/fatal. |
| **modulos_f90/run2** | Añadir PBC_Mod.f90 y GeomUtils.f90 al inicio de la lista de compilación. |
| **modulos_f90/Main.f90** | Sustituir llamada a print_params() por print_summary() (o alias). |

---

## 4. Formato de input.txt esperado (después del cambio)

Líneas de la forma:

```
clave valor
```

o

```
clave valor1 valor2 ...
```

- Comentarios con `!`; líneas vacías ignoradas.
- Claves (sin importar mayúsculas/minúsculas) según set_int/set_real/set_char (p, dp, sigmetano, eps, acelx, acely, acelz, angux, anguy, anguz, acel, diel, t, bcx, bcy, bcz, dim, mat, nam, isot, ijpasos, ikpasos, mult2, ensemble, ncellmat, nestado, canonicalmolecules).
- Se debe proporcionar un `input.txt` de ejemplo en el nuevo formato y, si se desea, una nota o script para convertir el input antiguo (orden fijo) al nuevo.

---

*Documento preparado para aplicar en Fase 1 (input clave–valor) y preparación Fase 2 (celda en InputParams). No se han modificado aún los fuentes.*
