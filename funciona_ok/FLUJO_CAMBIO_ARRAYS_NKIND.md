# Flujo de trabajo — cambio de arrays por `NKIND`

Este documento describe **el flujo de trabajo**, **los casos donde se rompe la lógica actual**, y **cómo se repararía** el problema de dimensionado cuando hay múltiples especies con tipos atómicos distintos.  
No se aplican cambios de código aquí.

---

## 1) Flujo de trabajo propuesto

1. **Identificar arrays dependientes de “tipos atómicos”**  
   Estos deben dimensionarse con `NKIND` (número total de tipos definidos en `LJ.dat`).

2. **Identificar arrays dependientes de “átomos por molécula”**  
   Estos deben dimensionarse con `maxAtoms` (máximo número de átomos en una especie).

3. **Leer `NKIND` antes de reservar arrays globales**  
   Se requiere leer la primera línea de `LJ.dat` antes de alocar `EPSI`, `SIGM`, `Q` y `USS`.

4. **Separar dimensiones y realocar**  
   - `EPSI`, `SIGM`, `Q` → tamaño `NKIND`  
   - `USS` → `USS(5000, NKIND, NKIND)`  
   - `NATOMKIND`, `NSYM`, `RX0/RY0/RZ0`, `RX/RY/RZ` → tamaño `maxAtoms`

5. **Agregar validaciones**  
   - Verificar que `max(NATOMKIND)` ≤ `NKIND`  
   - Verificar que todos los tipos leídos estén en `1..NKIND`

6. **Probar con casos representativos**  
   - Una sola especie  
   - Mezcla con tipos disjuntos  
   - Mezcla con tipos compartidos  

---

## 2) Casos donde hoy se rompe la lógica

### Caso A: Mezcla con tipos disjuntos
- Molécula 1: **A‑B** (2 átomos, tipos A y B)  
- Molécula 2: **C‑D‑F** (3 átomos, tipos C, D, F)  
- `maxAtoms = 3`  
- `NKIND = 5`

**Qué ocurre hoy**  
`EPSI/SIGM/Q/USS` se dimensionan con `maxAtoms = 3`, pero se accede con índices `1..5`.  
Resultado: acceso fuera de rango o memoria corrupta.

---

### Caso B: Varias especies pequeñas con tipos distintos
- 4 especies, cada una con 2 átomos  
- Tipos atómicos distintos totales: 6  
- `maxAtoms = 2`, `NKIND = 6`

**Qué ocurre hoy**  
Los arrays por tipo quedan en tamaño 2, pero se usan índices hasta 6.

---

### Caso C: Tipo atómico fuera de rango
- `NATOMKIND` contiene un tipo `8`  
- `LJ.dat` define `NKIND = 6`

**Qué ocurre hoy**  
Se accede a `EPSI(8)` y `USS(:,8,*)` que no existen.

---

## 3) Cómo se reparan esos casos

### Reparación base (cambio de dimensionado)
1. **Separar dimensiones**  
   - Arrays de tipos → `NKIND`  
   - Arrays de átomos por molécula → `maxAtoms`

2. **Mover lectura de `NKIND`**  
   - Leer `LJ.dat` antes de alocar `EPSI/SIGM/Q/USS`.

3. **Alocar correctamente**  
   - `EPSI(NKIND)`, `SIGM(NKIND)`, `Q(NKIND)`  
   - `USS(5000, NKIND, NKIND)`

---

### Validaciones mínimas
- Verificar que **todos** los tipos `NATOMKIND` están en `1..NKIND`.  
  Si no, emitir error y detener.

---

## 4) Resumen

- El problema aparece cuando **`NKIND > maxAtoms`**, algo común en mezclas con tipos disjuntos.  
- La reparación correcta es **dimensionar por `NKIND`** los arrays de tipos y dejar `maxAtoms` para los arrays por molécula.  
- Con validaciones simples, se evita que entradas inconsistentes produzcan fallos silenciosos.

