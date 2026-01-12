# MEJORAS Y REDUNDANCIAS EN Main.f90

## 🔴 PROBLEMAS CRÍTICOS

### 1. **Línea 433: Parámetros incorrectos en `call out()`** ⚠️ ERROR DE COMPILACIÓN
```fortran
call out(temp, z, sigma, eps, rcut, v, va, vg, w, ghost, jpasos, &
         nmin, nmaxi, canonicalmolecules)
```
**Problema CRÍTICO**: 
- La subrutina `OUT` NO acepta `nmin` y `nmaxi` como parámetros (ver `Out.f90` línea 20)
- La firma correcta es: `OUT(TEMP, Z, SIGMA, EPS, RCUT, V, VA, VG, W, GHOST, JPASOS, CANONICALMOLECULES)`
- Esto causará un **error de compilación**

**Solución**: Eliminar `nmin, nmaxi` de la llamada:
```fortran
call out(temp, z, sigma, eps, rcut, v, va, vg, w, ghost, jpasos, &
         canonicalmolecules)
```

### 2. **Línea 346: Índice fuera de bucle**
```fortran
anprom(I) = 0
ntotal(i) = 0
```
**Problema**: `I` está fuera del bucle `do i = 1, nmolec2` (que termina en línea 332).  
**Solución**: Debería ser un bucle explícito o usar el último valor de `I`:
```fortran
do I = 1, NMOLEC
   anprom(I) = 0
   ntotal(I) = 0
end do
```

### 3. **Líneas 623-624: Cierre de archivos no abiertos**
```fortran
close(21)
close(22)
```
**Problema**: Los archivos 21 y 22 nunca se abren en Main.  
**Solución**: Eliminar estas líneas o verificar si deberían abrirse.

---

## ⚠️ VARIABLES DECLARADAS PERO NO UTILIZADAS

> Nota: Muchas variables globales viven en los módulos (`SimulationData`, `AdsorbateInput`, `PhysicalConstants`). Las siguientes son locales en `Main` y no se usan allí. Revisa antes de borrarlas si piensas mover lógica aquí.

### Variables locales no usadas:
- **Línea 87**: `NMAX` - nunca se usa
- **Línea 89**: `RXNEW, RYNEW, RZNEW` - nunca se usan (solo se pasan a subrutinas)
- **Línea 92-93**: `DELTV, DELTW` - solo se usan localmente en un bloque
- **Línea 97**: `RMIN` - nunca se usa
- **Línea 99**: `OVRLAP` - nunca se usa
- **Línea 102**: `MOLEC1` - nunca se usa
- **Línea 113**: `NMATOM` - nunca se usa
- **Línea 115**: `IKIND, NS` - nunca se usan
- **Línea 141**: `NCONFMIN, NCONFMAX` - nunca se usan

**Recomendación**: Eliminar estas declaraciones o comentarlas si se planean usar en el futuro.

---

## 🔧 REDUNDANCIAS Y MEJORAS

### 1. **Línea 197: Asignación redundante**
```fortran
P = P ! * 1333.22
```
**Problema**: Asignación innecesaria con comentario.  
**Solución**: Eliminar la línea o descomentar si es necesario:
```fortran
! P = P * 1333.22  ! Conversión si es necesaria
```

### 2. **Líneas 404 y 416: Asignación duplicada de MULT**
```fortran
MULT = 1
...
MULT = 1
```
**Problema**: `MULT` se asigna dos veces.  
**Solución**: Eliminar la primera asignación (línea 404) ya que se sobrescribe en 416:
```fortran
! Eliminar línea 404, mantener solo:
MULT = 1
if (JPASOS.eq.1) then
   MULT = mult2
end if
```

### 3. **Línea 563: Asignación redundante dentro del bucle**
```fortran
do I = 1, NC
   read(49,*) RXAI, RYAI, RZAI, EPSAI, SGCI, QACI, SYMBOL2
   write(40,*) SYMBOL2, RXAI, RYAI, RZAI
   ntotalGRAF = NC  ! ← Redundante, ya se asignó en línea 545
end do
```
**Problema**: `ntotalGRAF = NC` se asigna en cada iteración innecesariamente.  
**Solución**: Eliminar del bucle (ya está asignado antes).

### 4. **Líneas 555-558: Variables escalares redundantes**
```fortran
escalax = acelx
escalay = acely
escalaz = acelz
write(*,*) escalax, escalay, escalaz
```
**Problema**: Se crean variables intermedias innecesarias.  
**Solución**: Usar directamente:
```fortran
write(*,*) acelx, acely, acelz
```

### 5. **Línea 533: Variable ESCALA redundante**
```fortran
ESCALA = acel
```
**Problema**: Se puede usar `acel` directamente.  
**Solución**: Reemplazar `ESCALA` por `acel` en línea 571-573:
```fortran
RXN = RX(jin,k,I)*acel
RYN = RY(jin,k,I)*acel
RZN = RZ(jin,k,I)*acel
```

### 6. **Líneas 583-585: Bloque vacío**
```fortran
if (ensemble2.eq.3) then
   ! (colocar aquí lógica adicional si aplica)
end if
```
**Problema**: Bloque condicional vacío.  
**Solución**: Eliminar si no se va a usar, o implementar la lógica.

---

## 📝 MEJORAS DE ESTILO Y ORGANIZACIÓN

### 1. **Líneas 148-150: Líneas en blanco excesivas**
```fortran
integer           auxmat



```
**Solución**: Eliminar líneas en blanco innecesarias.

### 2. **Línea 187: Formato de bucle obsoleto**
```fortran
do  i = 2, isot
```
**Solución**: Usar formato moderno:
```fortran
do i = 2, isot
```

### 3. **Línea 424: Uso de `goto` con etiquetas numéricas**
```fortran
IJ = ranf(DUMMY)*3 + 1
if (ensemble.eq.0) goto 30
goto (10,20,30,35) IJ
```
**Problema**: Uso de `goto` y etiquetas numéricas (estilo antiguo).  
**Solución**: Considerar usar `select case` o `if-else if`:
```fortran
IJ = int(ranf(DUMMY)*4) + 1
if (ensemble == 0) then
   call move(...)
else
   select case (IJ)
   case (1)
      call in(...)
   case (2)
      call out(...)
   case (3)
      call move(...)
   case (4)
      call change(...)
   end select
end if
```

### 4. **Línea 402: Unidad de archivo no declarada**
```fortran
write(58,*) aitest76
```
**Problema**: La unidad 58 no se abre explícitamente.  
**Solución**: Verificar si debería abrirse o eliminar si es debug.

### 5. **Línea 410: Apertura de archivo dentro del bucle JPASOS**
```fortran
open(unit=51, file='Ener'//CONFIG//'.TXT')
```
**Problema**: Se abre el archivo en cada iteración de JPASOS pero se cierra al final.  
**Solución**: Mover la apertura fuera del bucle JPASOS (antes de línea 400) y cerrar al final.

---

## 🗑️ DEALLOCATE: Verificaciones necesarias

### 1. **Línea 635: `deallocate(qac)`**
**Verificación**: `QAC` está en `SimulationData` como `allocatable`, así que está bien, pero verificar que se haya alocado.

### 2. **Línea 644: `deallocate(NATOM)`**
**Problema**: `NATOM` viene de `AdsorbateInput`, no de `SimulationData`.  
**Solución**: Verificar si debe desalocarse aquí o en el módulo.

### 3. **Líneas 641-644: Desalocación de arrays de módulos**
```fortran
deallocate(Z)
deallocate(X)
deallocate(N)
deallocate(NATOM)
```
**Problema**: `X`, `NATOM` vienen de `AdsorbateInput`, `N` de `SimulationData`.  
**Solución**: Verificar la propiedad de estos arrays y desalocarlos en el módulo correspondiente o aquí según diseño.

---

## ✅ RESUMEN DE ACCIONES RECOMENDADAS

### Críticas (deben corregirse):
1. ✅ Corregir `nmin, nmaxi` → `NMIN2, NMAXI2` en línea 433
2. ✅ Corregir índice `I` fuera de bucle en línea 346
3. ✅ Eliminar `close(21)` y `close(22)` o verificar apertura

### Importantes (mejoran el código):
4. ✅ Eliminar variables no utilizadas
5. ✅ Eliminar asignaciones redundantes (MULT, ntotalGRAF, escalax/y/z)
6. ✅ Mover apertura de archivo 51 fuera del bucle JPASOS
7. ✅ Reemplazar `goto` por `select case` o `if-else`

### Opcionales (mejoras de estilo):
8. ✅ Eliminar líneas en blanco excesivas
9. ✅ Modernizar formato de bucle `do`
10. ✅ Eliminar bloque condicional vacío (ensemble2)

---

## 📊 IMPACTO ESTIMADO

- **Reducción de líneas**: ~15-20 líneas eliminadas
- **Mejora de legibilidad**: Alta
- **Riesgo de cambios**: Bajo (solo correcciones críticas)
- **Mejora de mantenibilidad**: Media-Alta

---

## ✅ ESTADO DE REVISION Y CAMBIOS REALIZADOS

### Revisadas y realizadas
1. ✅ Corregir `nmin, nmaxi` en `call out(...)` (aplicado en `Main.f90`).
2. ✅ Corregir índice `I` fuera de bucle en línea 346 (aplicado en `Main.f90`).
3. ✅ Eliminar `close(21)` y `close(22)` (aplicado en `Main.f90`).

### Revisadas y no aplicadas
- (pendiente)
