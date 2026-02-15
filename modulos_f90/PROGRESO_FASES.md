# Progreso de Fases — Plan PBC/Input en modulos_f90

Documento de seguimiento del avance del plan definido en `PLAN_TRABAJO_PBC_INPUT.md`.
Se actualiza a medida que se completan o avanzan las fases.

---

## Resumen de estado

| Fase | Descripción | Estado | Notas |
|------|------------|--------|-------|
| Pre-fase | Revisión de PBC_Mod y GeomUtils | **Completada** | Bug min_image corregido; informe en `INFORME_PREFASE_PBC_GEOMUTILS.md` |
| 0 | Preparación y criterios | **Completada** | Respaldo en `respaldo/`; criterios definidos en el plan |
| 1 | Input clave-valor | **Completada** | `InputParams.f90` reescrito con parser clave-valor, `set_int/set_real/set_char`, `lowercase` |
| 2 | PBC_Mod + GeomUtils + Cell | **Completada** | Módulos integrados; `cell` y `cellR` construidos en InputParams; `run2` actualizado |
| 3 | Ángulos y dim en input | **Completada** | `angux/anguy/anguz`, `dim_flag`; `cell_from_lengths_angles` en InputParams |
| 4 | Potenciales usan Cell | **Completada** | `Potin.f90` y `Potout.f90` usan `min_image(cellR, ...)` |
| 5 | Movimientos MC usan Cell | **Completada** | `In.f90`, `Move.f90`, `change.f90` usan `cart_to_frac/wrap_by_pbc/frac_to_cart` con `cellR` |
| 6 | Estructura, estadística, salidas | **Completada** | cell_volume, fraccionales en Estructura/estadistica, VOL y XYZ en Main |
| 7 | Limpieza y triclínico | **Completada** | Variables muertas eliminadas, checks en fraccionales, tests triclínicos |

---

## Detalle por fase

### Pre-fase — Completada
- Revisión de `PBC_Mod.f90` y `GeomUtils.f90` de `roto/`.
- Bug detectado y corregido en `min_image`: usaba `cart_to_frac` con `centered=.true.` que envolvía ejes no periódicos. Corrección: `matmul(c%Ainv, r2-r1)` + `where (c%pbc) s = s - nint(s)`.
- Informe completo en `INFORME_PREFASE_PBC_GEOMUTILS.md`.

### Fase 0 — Completada
- Respaldo del código original en directorio `respaldo/`.
- Criterios de éxito definidos: compilación correcta + resultados numéricos equivalentes.

### Fase 1 — Completada
- `InputParams.f90`: parser clave-valor con `read_input('input.txt')`.
- Rutinas internas: `set_int`, `set_real`, `set_char`, `lowercase`, `info`, `warn`, `fatal`.
- `input_ejemplo.txt` disponible como referencia del formato nuevo.
- Todas las variables públicas mantienen los mismos nombres para compatibilidad.

### Fase 2 — Completada
- `PBC_Mod.f90` y `GeomUtils.f90` copiados/adaptados e integrados en `modulos_f90/`.
- `InputParams.f90` construye `cell` (Angstroms) y `cellR` (unidades reducidas, dividido por ACEL).
- `run2` actualizado: compila PBC_Mod → GeomUtils → InputParams → resto.
- Test: `test_pbc.f90` verifica celda ortorrómbica, min_image, wrap_by_pbc, celda triclínica.

### Fase 3 — Completada
- Claves `angux/alpha`, `anguy/beta`, `anguz/gamma`, `dim` disponibles en input.
- Valores por defecto: 90° y dim=3 (ortorrómbico por defecto).
- `cell_from_lengths_angles` construye la celda con ángulos y dimensionalidad.

### Fase 4 — Completada
- **Potin.f90**: Usa `min_image(cellR, r1, r2)` para calcular vectores mínima imagen en el cálculo de energía adsorbato-adsorbato al insertar.
- **Potout.f90**: Mismo esquema para el cálculo al eliminar.
- Eliminadas las fórmulas inline `RXIJ - BCX*XMAX*ANINT(RXIJ/XMAX)`.

### Fase 5 — Completada
- **In.f90**: Usa `cart_to_frac(cellR, pos)` → `wrap_by_pbc(...)` → `frac_to_cart(cellR, s)` para wrap de posiciones de prueba.
- **Move.f90**: Mismo patrón para rotación (CASE 1) y traslación (CASE 2).
- **change.f90**: Mismo patrón.
- Patrón establecido: todas las conversiones pasan por `cellR` (celda en unidades reducidas).

### Fase 6 — Completada
**Archivos modificados:** `PBC_Mod.f90`, `Estructura.f90`, `estadistica.f90`, `Main.f90`.

**Cambios realizados:**
- **PBC_Mod.f90**: Añadida función pública `cell_volume(c)` que calcula `det(A)` por cofactores.
- **Estructura.f90**: Importa `cell` de InputParams. Check "dentro de la caja" reemplazado: de `RXA >= -ACELX/2` a `all(abs(matmul(cell%Ainv, pos)) <= 0.5)` (fraccionales raw, sin wrap).
- **estadistica.f90**: Importa `cellR` y `cart_to_frac`. Índices CNF calculados vía `cart_to_frac(cellR, pos)` con clamp de rango para evitar out-of-bounds.
- **Main.f90**: `VOL = cell_volume(cellR)`. Salida XYZ usa `frac_to_cart(cell, cart_to_frac(cellR, pos))`. XMAX/YMAX/ZMAX mantenidos temporalmente para impresión.

**Verificación**: Compilación correcta con `bash run2` (solo warnings pre-existentes). Para ortorrómbico, todas las conversiones son algebraicamente idénticas al código anterior.

### Fase 7 — Completada
**Archivos modificados:** `Main.f90`, `In.f90`, `Move.f90`, `change.f90`, `test_pbc.f90`, `input_ejemplo.txt`.

**Cambios realizados:**
- **Main.f90**: Eliminadas variables muertas `XMAX/YMAX/ZMAX` (declaración y asignación).
- **In.f90**: Generación de posiciones aleatorias via fraccionales uniformes `[-0.5, 0.5)` + `frac_to_cart(cellR, s)`. Eliminadas `XMAX/YMAX/ZMAX`. Check inside-box en fraccionales (solo ejes no periódicos).
- **Move.f90**: Checks inside-box en fraccionales (2 bloques: rotación y traslación).
- **change.f90**: Check inside-box en fraccionales (1 bloque).
- **test_pbc.f90**: Tests 7-9 añadidos — min_image triclínico (|dr|²=s^T·G·s), cell_volume vs fórmula analítica, round-trip reducidas→Å.
- **input_ejemplo.txt**: Ejemplo comentado de celda triclínica.

**Verificación**: Compilación OK con `bash run2`. 10/10 tests pasan (incluyendo 3 nuevos triclínicos). Para ortorrómbico, comportamiento idéntico: `frac_to_cart(cellR, s)` con s∈[-0.5,0.5) produce las mismas posiciones que `s*XMAX`, y checks en fraccionales son equivalentes a `abs(RX1/XMAX) > 0.5`.

---

## Sistemas de coordenadas en el código

| Variable | Archivo | Unidades | Rango típico | Propósito |
|----------|---------|----------|--------------|-----------|
| RXA, RYA, RZA | Estructura.f90 | Angstroms | [-ACEL/2, ACEL/2] | Coords. crudas de superficie |
| RXC, RYC, RZC | SimulationData | Angstroms escalados (×σ/σ_CH4) | Similar | Coords. de superficie para potencial |
| RX, RY, RZ | SimulationData | Unidades reducidas (cart. en cellR) | [-0.5, 0.5] | Posiciones adsorbatos (almacenamiento principal) |
| RX1, RY1, RZ1 | SimulationData | Unidades reducidas | [-0.5, 0.5] | Coords. temporales en movimientos MC |
| initconf.txt | Main.f90 | Unidades reducidas | [-0.5, 0.5] | Archivo de reinicio |
| CONFIG*.xyz | Main.f90 | Angstroms (RX×ACEL) | Espacio real | Salida XYZ para visualización |

**Relaciones clave:**
- Angstroms → Reducidas: `r_red = r_Å / ACEL`
- Reducidas → Fraccionales de cellR: `s = cart_to_frac(cellR, r_red)`
- Fraccionales → Angstroms: `r_Å = frac_to_cart(cell, s)`
- Volumen reducido: `VOL = det(cellR%A)` (para ortorrómbico = XMAX×YMAX×ZMAX)

---

*Última actualización: febrero 2026. Referencia: PLAN_TRABAJO_PBC_INPUT.md*
