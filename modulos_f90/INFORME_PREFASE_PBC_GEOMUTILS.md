# Informe Pre-fase: revisión de PBC_Mod y GeomUtils

Revisión de los módulos **PBC_Mod.f90** y **GeomUtils.f90** (versión en `roto/`) según el plan PLAN_TRABAJO_PBC_INPUT.md, para validar lógica, coherencia de algoritmos y ausencia de bugs antes de integrarlos en modulos_f90.

---

## 1. PBC_Mod.f90

### 1.1 Convención de vectores (cell_from_lengths_angles)

- **Convención**: a lo largo de x, b en el plano xy, c con ángulos α (b∧c), β (a∧c), γ (a∧b).
- **Fórmulas**:
  - `ax = [a, 0, 0]`
  - `bx = [b*cos(γ), b*sin(γ), 0]`
  - `cx`: componente x = c·cos(β); componente y = c·(cos(α)−cos(β)cos(γ))/sin(γ); componente z por normalización (raíz).
- **Ortorrómbico** (α=β=γ=90°): cosa=cosb=cosg=0, sing=1 → cx = (0, 0, c). No hay división por cero.
- **Conclusión**: Coherente con convención cristalográfica estándar.

### 1.2 Inversa inv3 (cofactores)

- Fórmula por cofactores: `invM(i,j) = C(j,i)/det`, con C = matriz de cofactores de M.
- Revisión elemento a elemento: signos e índices correctos.
- Determinante: `det = M(1,1)*invM(1,1) + M(1,2)*invM(2,1) + M(1,3)*invM(3,1)` (desarrollo por fila 1).
- **Nota**: No hay comprobación de det≠0; si la celda es degenerada (ángulos 0°/180° o longitudes cero) puede haber fallo. Queda documentado; opcional añadir guarda o comprobación en fases posteriores.
- **Conclusión**: Fórmula correcta.

### 1.3 cart_to_frac / frac_to_cart

- **cart_to_frac**: s = Ainv·r; si `centered=.true.` se hace s = s − nint(s) → fraccionales en [-0.5, 0.5). Coherente.
- **frac_to_cart**: r = A·s. Correcto.
- **Conclusión**: Coherente con el uso en el resto del código.

### 1.4 wrap_frac_centered y wrap_cart

- **wrap_frac_centered**: sw = s − nint(s) para los tres ejes; no usa `pbc` (envuelve siempre los tres). Documentado en el plan.
- **wrap_cart**: cadena cart → frac → wrap_frac_centered → frac_to_cart; no usa `pbc`. Para wrap respetando PBC debe usarse otro flujo (p. ej. wrap_by_pbc en GeomUtils + frac_to_cart).
- **Conclusión**: Comportamiento correcto y documentado.

### 1.5 min_image(c, r1, r2) — bug corregido

- **Problema**: Se usaba `s = cart_to_frac(c, r2 − r1)`. Con `centered=.true.`, `cart_to_frac` aplica s = s − nint(s) en **los tres** ejes. Para geometría tipo slab (p. ej. pbc=(T,T,F)) el eje z no es periódico y el vector mínima-imagen debe conservar la diferencia real en z; al aplicar nint(s3) se envolvía incorrectamente el eje no periódico.
- **Corrección aplicada** (en `roto/PBC_Mod.f90`): Calcular el vector diferencia en fraccionales sin usar `cart_to_frac`, y envolver solo donde `c%pbc`:
  - `s = matmul(c%Ainv, r2 - r1)`
  - `where (c%pbc) s = s - nint(s)`
  - `dr_cart = frac_to_cart(c, s)`
- **Conclusión**: Bug corregido; min_image respeta correctamente ejes no periódicos.

### 1.6 pbc_ORTO(dr, lbox)

- Fórmula: `dr_min = dr − lbox*nint(dr/lbox)`. Válida para caja ortorrómbica; lbox = [Lx, Ly, Lz] en el mismo orden que el resto de rutinas.
- **Conclusión**: Correcto.

### 1.7 Kind rk

- `selected_real_kind(12, 99)` usado de forma consistente en el módulo.
- **Conclusión**: OK.

---

## 2. GeomUtils.f90

### 2.1 cell_to_metric

- **G = transpose(A)·A**: métrica G_ij = Σ_k A(k,i)A(k,j); G simétrica; r² = s^T G s para vector fraccional s.
- Devuelve Ainv y los tres flags pbc (pbcx, pbcy, pbcz) desde la celda.
- **Conclusión**: Correcto; G es la misma que s^T G s en r2_min_image_frac.

### 2.2 r2_min_image_frac(G, s1, s2, s3)

- Fórmula: r² = G(1,1)*s1² + G(2,2)*s2² + G(3,3)*s3² + 2*(G(1,2)*s1*s2 + G(1,3)*s1*s3 + G(2,3)*s2*s3). Coincide con s^T G s para s = (s1,s2,s3). Índices correctos; G es transpose(A)*A.
- **Conclusión**: Correcto.

### 2.3 wrap_by_pbc(s1, s2, s3, pbcx, pbcy, pbcz)

- Aplica s = s − nint(s) solo en los ejes donde el flag es true. Para slab (pbcz=.false.) s3 no se modifica. Argumentos inout (s1, s2, s3) se actualizan correctamente.
- **Conclusión**: Crítico para slab; comportamiento correcto.

### 2.4 Dependencia del tipo Cell

- GeomUtils usa `Cell_t => Cell` desde PBC_Mod. Compatible con PBC_Mod.
- **Conclusión**: OK.

---

## 3. Checklist de posibles bugs (plan)

| Ítem | Estado |
|------|--------|
| Celda ortorrómbica (α=β=γ=90°): A diagonal en bloque, Ainv y G coherentes | Verificado; sin división por cero |
| Determinante de A: no nulo en casos típicos; manejo/documentación de degenerados | Documentado; sin guarda explícita |
| min_image con pbc(3)=[T,T,F]: z no se envuelve; x,y sí; resultado en cartesianas correcto | Corregido (cálculo directo + where (c%pbc)) |
| Consistencia ortorrómbica: min_image vs dr−L*nint(dr/L) en ejes periódicos | Coherente con la corrección |
| wrap_by_pbc: ejes pbc=.false. no cambian s; .true. → [-0.5, 0.5) | Verificado |
| Uso de nint para caja centrada s−nint(s) | Correcto en todos los usos |

---

## 4. Resumen y criterio de cierre

- **Revisiones realizadas**: Convención de vectores, inv3, cart_to_frac/frac_to_cart, wrap_frac_centered, wrap_cart, min_image, pbc_ORTO, cell_to_metric, r2_min_image_frac, wrap_by_pbc.
- **Bug detectado y corregido**: En **min_image**, uso de `cart_to_frac` para el vector diferencia con `centered=.true.` envolvía también ejes no periódicos; corregido en `roto/PBC_Mod.f90` usando `matmul(c%Ainv, r2-r1)` y `where (c%pbc) s = s - nint(s)`.
- **Entregables**: Este informe; versión corregida de PBC_Mod en `roto/` lista para reutilizar al integrar en modulos_f90 (Fase 2).

**Criterio de cierre Pre-fase**: Revisión hecha; bug de min_image resuelto. Se puede pasar a Fase 1 (input clave–valor) y luego Fase 2 (integrar PBC_Mod y GeomUtils en modulos_f90 usando la versión de `roto/` actualizada).

---

*Fecha: febrero 2025. Referencia: PLAN_TRABAJO_PBC_INPUT.md, Pre-fase.*
