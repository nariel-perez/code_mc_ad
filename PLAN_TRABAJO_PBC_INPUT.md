# Plan de trabajo: input tipo roto y PBC extendidas en modulos_f90

Objetivo: incorporar en **modulos_f90** (1) el formato de entrada clave–valor usado en **roto** y (2) la lógica de PBC modular y extensible (ortorrómbico → triclínico), avanzando por fases sin romper lo que ya funciona.

**Referencia**: diseño actual en `roto/` (InputParams, PBC_Mod, GeomUtils, uso en Potin/Potout/In/Move, etc.).  
**Base**: código actual en `modulos_f90/` (input por orden fijo, PBC con `BCX/BCY/BCZ` y fórmulas inline). Este código está documentado (DIAGRAMA_FLUJO.md, ESTADO_DEL_CODIGO.md) y se tomarán **copias de respaldo** antes de reestructurar, de modo que siempre se pueda volver a una versión estable.

---

## Pre-fase: Análisis y validación de PBC_Mod y GeomUtils

**Objetivo**: Revisar con cuidado los módulos **PBC_Mod.f90** y **GeomUtils.f90** (en `roto/` o en la versión que se vaya a reutilizar) para detectar y corregir posibles bugs o artefactos antes de integrarlos en modulos_f90. Evitar que errores en la geometría o en el wrap contaminen la simulación.

### Contenido del análisis

- **PBC_Mod.f90**
  - **Convención de vectores**: verificar que la matriz `A` (columnas = vectores a, b, c) y la construcción en `cell_from_lengths_angles` siguen la convención cristalográfica estándar (a a lo largo de x, b en el plano xy, c con ángulos α, β, γ). Revisar fórmulas de ax, bx, cx (en particular el término de c que involucra raíz cuadrada) para ángulos no triviales.
  - **Inversa inv3**: revisar la fórmula por cofactores (signos, índices) y que el determinante se use correctamente (det ≠ 0). Casos degenerados (ángulos 0° o 180°, longitudes cero) pueden dar det = 0; documentar o añadir guardas si aplica.
  - **cart_to_frac / frac_to_cart**: s = Ainv·r y r = A·s; comprobar que con `centered=.true.` el ajuste s = s - nint(s) en cart_to_frac lleva a fraccionales en [-0.5, 0.5). Coherencia con el uso en el resto del código (roto usa centered en la construcción).
  - **wrap_frac_centered**: s - nint(s) para los tres ejes; correcto para “caja centrada”. Dejar claro que **no** respeta `pbc` (envuelve siempre los tres ejes).
  - **wrap_cart**: cadena cart → frac → wrap_frac_centered → frac_to_cart; no usa `pbc`. Documentar que para “wrap respetando PBC” se debe usar otro flujo (p. ej. wrap_by_pbc en fraccionales y luego frac_to_cart).
  - **min_image(c, r1, r2)**: comprobar que `where (c%pbc) s = s - nint(s)` se aplica solo a los ejes periódicos y que en ejes no periódicos no se modifica s (vector mínimo imagen correcto para slab 1 1 0). Verificar que el resultado en cartesianas (frac_to_cart(c, s)) es consistente.
  - **pbc_ORTO(dr, lbox)**: dr_min = dr - lbox*nint(dr/lbox); solo válido para caja ortorrómbica; comprobar que lbox = [Lx, Ly, Lz] en el mismo orden que usan el resto de rutinas.
  - **Kind rk**: uso consistente de `selected_real_kind(12, 99)` para evitar mezclas de precisión con el resto del código.

- **GeomUtils.f90**
  - **cell_to_metric**: G = transpose(A)*A (métrica); comprobar que G es simétrica y que r² = s^T G s para un vector fraccional s. Devolver Ainv y los tres flags pbc desde la celda; revisar que no se copien por error solo dos ejes.
  - **r2_min_image_frac(G, s1, s2, s3)**: r² = G(1,1)*s1² + G(2,2)*s2² + G(3,3)*s3² + 2*(G(1,2)*s1*s2 + G(1,3)*s1*s3 + G(2,3)*s2*s3). Coincide con s^T G s para s = [s1,s2,s3]. Verificar índices y que G sea la misma que devuelve cell_to_metric (transpose(A)*A, no A*transpose(A)).
  - **wrap_by_pbc(s1, s2, s3, pbcx, pbcy, pbcz)**: aplicar s = s - nint(s) solo en los ejes donde el flag es true. Crítico para slab: si pbcz=.false., s3 no debe modificarse. Revisar que los argumentos inout (s1, s2, s3) se actualicen correctamente en todos los casos.

### Checklist de posibles bugs y artefactos

- [ ] Celda ortorrómbica (α=β=γ=90°): A es diagonal en bloque (a,0,0), (b,0,0), (c_x, c_y, c_z) con c_z ≠ 0; Ainv y G coherentes.
- [ ] Determinante de A: no nulo para ángulos y longitudes típicos; manejo (o documentación) de casos degenerados.
- [ ] min_image con pbc(3)=[T,T,F]: vector en z no se envuelve; en x e y sí; resultado en cartesianas correcto.
- [ ] Consistencia numérica: para ortorrómbico, min_image(c, r1, r2) debe coincidir con dr - L*nint(dr/L) componente a componente en los ejes periódicos (tolerancia por redondeo).
- [ ] wrap_by_pbc: ejes con pbc=.false. no deben cambiar el valor de s; ejes con .true. deben quedar en [-0.5, 0.5).
- [ ] Uso de nint vs floor/ceiling: en Fortran nint redondea al entero más cercano; para wrap a caja centrada s - nint(s) es lo correcto. Verificar que no haya desajustes en los límites ±0.5.

### Entregables Pre-fase

- **Informe breve** (o anotaciones en el plan / en un .md) con: (1) lista de revisiones hechas, (2) bugs o dudas encontradas y su resolución (o ticket para resolver antes de Fase 2), (3) confirmación de que inv3, min_image, wrap_by_pbc y r2_min_image_frac están validados para ortorrómbico y, si se probó, para un caso triclínico simple.
- **Opcional**: pequeño programa o bloque de tests (p. ej. en un directorio tests/ o al final de un módulo de prueba) que compruebe: celda ortorrómbica Lx,Ly,Lz vs min_image y vs pbc_ORTO; r² con r2_min_image_frac vs norma al cuadrado de min_image; wrap_by_pbc con distintos (pbcx,pbcy,pbcz). No es obligatorio para cerrar la pre-fase, pero refuerza la confianza.
- **Criterio de cierre**: No pasar a Fase 1 con la integración de estos módulos hasta tener la Pre-fase cerrada (revisión hecha y, si hubo cambios, versión “limpia” de PBC_Mod y GeomUtils lista para copiar a modulos_f90).

---

## Fase 0: Preparación y criterios

- **Respaldos**: Antes de aplicar cualquier cambio estructural, hacer (o verificar) **copias de respaldo** del directorio modulos_f90. Así siempre se puede volver a la versión estable y documentada (DIAGRAMA_FLUJO.md, ESTADO_DEL_CODIGO.md) si algo sale mal.
- **Criterio de éxito por fase**: compilación correcta y resultados numéricos equivalentes (mismo sistema, mismo input conceptual) antes y después del cambio de la fase, cuando aplique.
- **Orden**: Pre-fase (revisar PBC_Mod y GeomUtils); luego Fase 1 (input clave–valor); después introducción del tipo Cell y PBC (Fases 2–3); por último sustitución progresiva de PBC inline (Fases 4–7).
- **Compatibilidad**: hasta que se decida lo contrario, mantener la posibilidad de leer el **input por orden fijo** (legacy) además del nuevo formato clave–valor, o documentar claramente el nuevo formato y migrar los `input.txt` de prueba.

---

## Fase 1: Nuevo formato de input (clave–valor) en modulos_f90

**Objetivo**: que modulos_f90 pueda leer el mismo estilo de input que roto (líneas `clave valor` o `clave valor1 valor2...`), sin cambiar aún la geometría ni las PBC.

### 1.1 Diseño

- Añadir (o reemplazar en un módulo) un **parser de líneas**:
  - Ignorar líneas vacías y comentarios (p. ej. `!`).
  - Por línea: leer un **token** (clave) y luego uno o más valores (entero, real o cadena según clave).
  - Mapear claves a variables ya existentes en InputParams (P, dp, sigmetano, eps, ACEL, acelx, acely, acelz, diel, BCX, BCY, BCZ, T, mat, nam, isot, ijpasos, ikpasos, mult2, ensemble, NCELLMAT, NESTADO, canonicalmolecules).
- Definir **nombres de claves** compatibles con roto donde tenga sentido (p, dp, sigmetano, eps, acelx, acely, acelz, diel, bcx, bcy, bcz, t, mat, nam, isot, ijpasos, ikpasos, mult2, ensemble, ncellmat, nestado, canonicalmolecules; y opcionalmente acel, alpha/beta/gamma para fases posteriores).
- Rutinas auxiliares: `lowercase(trim(key))`, `set_int(key, val)`, `set_real(key, val)`, `set_char(key, val)` (o equivalente) que actualicen las variables del módulo según la clave.

### 1.2 Integración sin romper legacy

- **Opción A**: Añadir `read_input_keyval(fname)` en InputParams; en `read_input(fname)` detectar si el archivo empieza por un número (primera columna numérica) y llamar a la lectura por orden fijo o a `read_input_keyval` según corresponda.
- **Opción B**: Reemplazar la lectura por una sola rutina clave–valor y documentar el nuevo formato; proporcionar un `input.txt` de ejemplo en el nuevo formato y, si se desea, un script o nota para convertir inputs viejos.

### 1.3 Entregables Fase 1

- Input leído por clave–valor con el mismo conjunto de parámetros que hoy.
- Mismo comportamiento de la simulación (mismos resultados para el mismo sistema) usando un input equivalente en el nuevo formato.
- Documentación breve del formato (ejemplo de `input.txt` y lista de claves en un .md o en comentarios del código).

---

## Fase 2: Introducir PBC_Mod y el tipo Cell (sin cambiar aún la física)

**Objetivo**: tener en modulos_f90 el módulo **PBC_Mod** y el tipo **Cell** como en roto, y construir una celda ortorrómbica a partir de los parámetros actuales (acelx, acely, acelz, BCX, BCY, BCZ), sin sustituir todavía las PBC en el resto del código.

### 2.1 Archivos nuevos (copiar/adaptar desde roto)

- **PBC_Mod.f90** (o nombre análogo):
  - Tipo `Cell` con: longitudes a,b,c; ángulos α,β,γ; matrices A(3,3), Ainv(3,3); lógicos pbc(3); entero dim; opcional centered.
  - `cell_from_lengths_angles(c, a, b, c, alpha, beta, gamma [, dim] [, centered])`.
  - Para ortorrómbico: llamar con alpha=beta=gamma=90.
  - `cell_update(c)` que calcule Ainv (inv3).
  - Funciones puras: `cart_to_frac`, `frac_to_cart`, `wrap_frac_centered`, `wrap_cart`, `min_image(c, r1, r2)` respetando `c%pbc`.
  - Opcional: `pbc_ORTO(dr, lbox)` para uso legacy ortorrómbico.
  - Kind `rk` (selected_real_kind) para consistencia.

- **GeomUtils.f90** (o nombre análogo):
  - `cell_to_metric(c, G, Ainv, pbcx, pbcy, pbcz)`.
  - `r2_min_image_frac(G, s1, s2, s3)`.
  - `wrap_by_pbc(s1, s2, s3, pbcx, pbcy, pbcz)`.

### 2.2 Integración en InputParams (modulos_f90)

- Añadir en InputParams (o donde se lean parámetros): después de leer acelx, acely, acelz y BCX, BCY, BCZ, **construir** un objeto `cell`:
  - Por ahora: `cell_from_lengths_angles(cell, acelx, acely, acelz, 90._rk, 90._rk, 90._rk, dim=3, centered=.true.)` (o usar ACEL como referencia de longitud según cómo esté definida la caja en Main).
  - Asignar `cell%pbc(1:3) = [ BCX/=0, BCY/=0, BCZ/=0 ]`, `cell%dim = 3` (o un dim_flag si ya se leyó en Fase 1).
- InputParams debe `use PBC_Mod` y exponer `cell` (y si aplica `rk`) para el resto del programa.
- **No** cambiar todavía In.f90, Move.f90, Potin.f90, Potout.f90, Potencial.f90, etc.: que sigan usando BCX, BCY, BCZ, XMAX, YMAX, ZMAX como hasta ahora.

### 2.3 Entregables Fase 2

- PBC_Mod y GeomUtils compilando e integrados en el build (run2 o equivalente).
- InputParams construye `cell` después de leer; el resto del código sigue igual y los resultados se mantienen.
- Breve verificación: en Main (o un test mínimo), comprobar que para la celda ortorrómbica actual, `min_image(cell, r1, r2)` coincide con el vector mínima-imagen que se obtendría con la fórmula actual (dr - L*nint(dr/L)) en cada eje periódico.

---

## Fase 3: Leer ángulos y dim en el input; celda desde longitudes y ángulos

**Objetivo**: que el input pueda definir una celda por longitudes **y** ángulos (α, β, γ), y opcionalmente dim (2/3), para preparar triclínico sin usarlo aún en potenciales/movimientos.

### 3.1 Input

- Añadir claves (si no existen): por ejemplo `acel`, `angux`/`alpha`, `anguy`/`beta`, `anguz`/`gamma`, `dim`.
- Valores por defecto: alpha=beta=gamma=90, dim=3, para no romper inputs que solo tengan longitudes.
- En read_input (o keyval): después de leer, llamar `cell_from_lengths_angles(cell, a_len, b_len, c_len, alpha_deg, beta_deg, gamma_deg, dim=dim_flag, centered=.true.)` y luego `cell%pbc = [BCX/=0, BCY/=0, BCZ/=0]`, `cell%dim = dim_flag`.

### 3.2 Compatibilidad

- Si no se leen ángulos ni dim, usar 90° y 3 como hasta ahora; la celda sigue siendo ortorrómbica y el comportamiento idéntico.

### 3.3 Entregables Fase 3

- Input con claves opcionales para ángulos y dim; celda construida con `cell_from_lengths_angles` en todos los casos.
- Tests rápidos: ortorrómbico (90,90,90) da los mismos A, Ainv que antes; variar α,β,γ y comprobar que A tiene la forma esperada (convención cristalográfica).

---

## Fase 4: Sustituir PBC inline por uso de Cell y GeomUtils (potenciales)

**Objetivo**: que las **tablas de potencial** (Potencial.f90, PotencialFF.f90 si aplica) y los cálculos de **energía** que usan distancias mínimas (Potin, Potout, Adpotin, Adpotout) pasen a usar la celda reducida `cellR` y las rutinas del módulo PBC/GeomUtils en lugar de fórmulas con BCX*XMAX*ANINT(...).

### 4.1 Estrategia

- En cada subrutina que hoy usa BCX, BCY, BCZ y XMAX, YMAX, ZMAX para imagen mínima:
  - Obtener una celda en unidades reducidas: `cellR = cell`, `cellR%A = cell%A / ACEL` (o el escalado que use Main: sigma, sigmetano, etc.), `cellR%update()`.
  - Llamar `cell_to_metric(cellR, G, Ainv, px, py, pz)`.
  - Para cada par de posiciones: Δr → Δs = Ainv * Δr (en fraccionales); llamar `wrap_by_pbc(s1,s2,s3, px,py,pz)`; calcular r² con `r2_min_image_frac(G, s1,s2,s3)`.
  - O, donde sea más simple, usar `min_image(cellR, r1, r2)` y luego norma al cuadrado.
- Eliminar (o dejar comentado temporalmente) las líneas que hacen `RXIJ = RXIJ - BCX*XMAX*ANINT(RXIJ/XMAX)` etc. en esos archivos.
- Asegurar que Potencial, Potin, Potout, Adpotin, Adpotout reciben o tienen acceso a `cell` (vía use InputParams o use PBC_Mod) y a ACEL (o el factor de escala correcto) para construir cellR.

### 4.2 Orden sugerido de archivos

1. **Potencial.f90** (adsorbato–superficie): sustituir cálculo de distancias por fraccionales + wrap_by_pbc + r2_min_image_frac (o min_image).
2. **Potin.f90** y **Potout.f90**: mismo esquema.
3. **Adpotin.f90** y **Adpotout.f90**: mismo esquema.
4. Si PotencialFF usa PBC, repetir el patrón ahí.

### 4.3 Entregables Fase 4

- Potenciales y energías calculados vía Cell + GeomUtils/PBC_Mod.
- Resultados numéricos de la simulación iguales (o dentro de tolerancia numérica) respecto a la versión anterior con PBC inline ortorrómbicas.

---

## Fase 5: Sustituir PBC en movimientos MC (In, Move, change)

**Objetivo**: que las **nuevas posiciones** y el **wrap** en In.f90, Move.f90 y change.f90 usen la celda y las utilidades (wrap_by_pbc, frac_to_cart, cart_to_frac) en lugar de fórmulas con BCX*XMAX*ANINT(...).

### 5.1 In.f90

- Tras generar posición de prueba y rotar: convertir a fraccionales con la celda reducida cellR; aplicar `wrap_by_pbc` solo en ejes con PBC; comprobar “inside” en ejes sin PBC (p. ej. que las fraccionales estén en [0,1) o [-0.5,0.5) según convención); volver a cartesianas con cellR%A * s.
- Reemplazar las líneas que hacen RX1(I) = RX1(I) - BCX*XMAX*ANINT(...) por este flujo (o por una rutina auxiliar que lo encapsule).

### 5.2 Move.f90 y change.f90

- Mismo criterio: posiciones trial en fraccionales → wrap_by_pbc → chequeo inside en ejes no periódicos → vuelta a cartesianas.
- Asegurar que el rechazo por “fuera de caja” en ejes sin PBC sea equivalente al actual (p. ej. RETURN si abs(coord) > 0.5 en ese eje).

### 5.3 Entregables Fase 5

- In, Move, change sin uso directo de BCX/BCY/BCZ en fórmulas de wrap; todo vía cell + GeomUtils/PBC_Mod.
- Misma aceptación/rechazo de movimientos para el mismo input ortorrómbico; resultados de la simulación equivalentes.

---

## Fase 6: Estructura y estadística; salidas (XYZ, CNF)

**Objetivo**: que **Estructura.f90** y **estadistica.f90** (y cualquier escritura que use posiciones o límites de caja) usen la celda y el wrap respetando PBC (wrap solo en ejes periódicos).

### 6.1 Estructura.f90

- Criterio de “dentro de la caja” y escritura del archivo truncado: usar cart_to_frac(cell, r), wrap solo donde cell%pbc, y condición de inside en ejes no periódicos (como en roto).
- Ajustar si hoy usa acelx, acely, acelz o ACEL directamente; pasar a usar cell (y si hace falta escalado en Å, usar las longitudes de cell en las unidades que correspondan).

### 6.2 estadistica.f90

- Asignación a celdas de la matriz CNF: posiciones en fraccionales, wrap_by_pbc; luego mapeo (s1,s2,s3) → índices (ICNF,JCNF,KCNF) de forma consistente con el rango [-NCELLMAT/2 : NCELLMAT/2] o el que se use.
- Asegurar que en ejes sin PBC solo se cuentan partículas dentro del primario (o documentar el criterio).

### 6.3 Main.f90: volumen, escritura XYZ, initconf

- Volumen: si hoy VOL = XMAX*YMAX*ZMAX, reemplazar por el volumen de la celda (det(A) o equivalente) en las unidades que use Main; para ortorrómbico debe coincidir con el valor actual.
- Escritura de CONFIG*.xyz y lectura del archivo truncado: seguir usando las posiciones en cartesianas; el archivo truncado ya vendrá filtrado por Estructura con el nuevo criterio.
- initconf: sin cambios de formato; solo asegurar que las posiciones que se escriben/leen son las mismas que usa el resto (cartesianas en la celda reducida o en Å según convenga).

### 6.4 Entregables Fase 6

- Estructura y estadística usando cell y wrap por PBC.
- Salidas (truncado, XYZ, CNF, initconf) coherentes con el nuevo esquema; resultados de simulación equivalentes para el caso ortorrómbico.

---

## Fase 7: Limpieza y opción triclínica

**Objetivo**: eliminar variables redundantes (si se puede dejar solo cell como fuente de verdad para geometría y PBC) y validar explícitamente un caso triclínico (α, β o γ ≠ 90°).

### 7.1 Limpieza

- Revisar que ningún módulo siga usando BCX, BCY, BCZ, XMAX, YMAX, ZMAX para PBC; si queda algo, migrarlo a cell/GeomUtils.
- Decidir si se mantienen acelx, acely, acelz (y ACEL) como variables de input que solo alimentan la construcción de cell, o si se leen solo longitudes/ángulos y cell es la única definición.
- Actualizar documentación (DIAGRAMA_FLUJO.md, ESTADO_DEL_CODIGO.md o similar) con el nuevo flujo de input y PBC.

### 7.2 Triclínico

- Añadir (o reutilizar) un input de prueba con α, β o γ ≠ 90° (por ejemplo un paralelogramo en 2D o una caja ligeramente inclinada).
- Comprobar que compila, que las distancias mínimas y el volumen son consistentes, y que la simulación corre sin errores; comparar si es posible con una referencia (analítica o otro código).

### 7.3 Entregables Fase 7

- Código sin PBC inline; toda la lógica de caja y periodicidad en PBC_Mod y GeomUtils.
- Documentación actualizada y al menos un caso triclínico verificado.
- Plan cerrado para “PBC extendidas” en modulos_f90; extensión futura (p. ej. más convenciones de celda) sobre esta base.

---

## Resumen del flujo de trabajo

| Etapa | Contenido principal | Criterio de éxito |
|-------|---------------------|-------------------|
| **Pre-fase** | Análisis y validación de PBC_Mod.f90 y GeomUtils.f90 | Revisión hecha; bugs/artefactos resueltos o documentados; listos para integrar |
| 0 | Preparación y criterios (input legacy/keyval, backups) | Criterios y orden definidos; respaldos de modulos_f90 hechos |
| 1 | Input clave–valor en modulos_f90 | Mismos parámetros leídos; misma simulación |
| 2 | PBC_Mod + GeomUtils + Cell ortorrómbico en InputParams | Compila; cell construido; resto sin cambios |
| 3 | Input con ángulos y dim; cell_from_lengths_angles | Ortorrómico por defecto; opción triclínico en input |
| 4 | Potenciales y energías usan cell + wrap_by_pbc + r2_min_image_frac | Mismos resultados que antes |
| 5 | In, Move, change usan cell y wrap/inside | Mismos resultados que antes |
| 6 | Estructura, estadística, volumen y salidas usan cell | Salidas coherentes; mismos resultados |
| 7 | Limpieza; validación triclínico; documentación | Sin PBC inline; triclínico verificado |

---

## Notas

- **Respaldos**: Antes de reestructurar, se harán copias de respaldo del código de modulos_f90 (bien documentado en DIAGRAMA_FLUJO.md y ESTADO_DEL_CODIGO.md) para poder volver atrás si hace falta.
- **Orden de compilación**: PBC_Mod primero (no depende de InputParams); GeomUtils depende de PBC_Mod; InputParams puede depender de PBC_Mod; el resto como hasta ahora. Actualizar run2 o el script de build para incluir los nuevos .f90.
- **Tests**: en cada fase, conservar al menos un input.txt y una corrida corta de referencia para comparar salidas (p. ej. número de partículas por isoterma, energías o calores en una o dos líneas de SALIDAACTIVADO).
- **Revertir**: si una fase introduce fallos, revertir solo los archivos de esa fase y mantener el resto; el plan está pensado para commits por fase.

Este documento es solo el **plan**; no se han modificado archivos del proyecto. Al implementar, seguir las fases en orden y validar cada una antes de pasar a la siguiente.
