MODULE global_var_automata
implicit none
! ====================================================================
! ====================================================================
!
! 				MODULO CON LAS VARIABLES GLOBALES
!					FUNCIONES Y SUBRUTINAS
!
! ====================================================================
! ====================================================================
! PARAMETROS
REAL(8), PARAMETER :: pi = 3.1415926535897932384626

! VARIABLES DE CARACTERES
CHARACTER(len=50) :: system_filename	! Fichero con el estado inicial
CHARACTER(len=800) :: cwd		! Directorio con los ficheros .f90
CHARACTER(len=800) :: save_dir 	! Directorio para salvar los resultados

! VARIABLES ENTERAS
INTEGER :: n_mat		! Tamanio del sistema (matriz n_mat x n_mat)
INTEGER :: n_phases		! Numero de fases
INTEGER :: masa			! Numero de celdas ocupadas por materia
INTEGER :: iter			! Iteracion actual
INTEGER :: iter_prom	! Numero de iteraciones para promedio energetico
INTEGER :: iter_min		! Numero de iteraciones minimas

! VARIABLES REALES
REAL(8) :: densidad		! Densidad para iniciar con conf aleatoria
REAL(8) :: T			! Temperatura
REAL(8) :: alpha_mas	! Coeficiente para pesar el aumento de masa
REAL(8) :: alpha_menos	! Coeficiente para pesar la disminucion de masa
REAL(8) :: crit_end		! Criterio de energia para estado final
REAL(8) :: E_inicial	! Energia incial (antes de evolucionar)
REAL(8) :: E_final		! Energia final (tras evolucionar)
REAL(8) :: E_1			! Energia 1 para el criterio de equilibrio
REAL(8) :: E_2			! Energia 2 para el criterio de equilibrio
REAL(8) :: delta_E_rel	! Variacion de energia relativa en la evolucion
REAL(8) :: mov			! Movilidad
REAL(8) :: porosidad	! Porosidad
REAL(8) :: E_minima 	! Energia superficial minima (toda la masa 
						! formando un circulo

! VARIABLES ARRAYS
REAL(8), DIMENSION(:,:), ALLOCATABLE :: kp		! Energias de enlace
												! por cada pareja de
												! fases
INTEGER, DIMENSION(:,:), ALLOCATABLE :: mat	! Sistema (matriz con 
												! los identificadores de
												! cada fase)
INTEGER, DIMENSION(:,:), ALLOCATABLE :: mat_prom! Sistema promedio
REAL(8), DIMENSION(:), ALLOCATABLE :: ener_prom	! Energias para promedio
												
! FUNCIONES Y SUBRUTINAS
CONTAINS

! *********************************************
! 				FUNCTION: KRON
! *********************************************
REAL(8) FUNCTION KRON(i_ind, j_ind)

! DECLARACION DE VARIABLES LOCALES
implicit none
INTEGER :: i_ind, j_ind		! Indices para la pareja de celdas

! FUNCION
kron = kp(i_ind, j_ind)

! FIN DE LA FUNCION
END FUNCTION KRON


! *********************************************
! 				FUNCTION: ENERGIA
! *********************************************
REAL(8) FUNCTION  ENERGIA(i_cell, j_cell, tipo)

! DECLARACION DE VARIABLES LOCALES
implicit none
INTEGER :: i_cell, j_cell	! Indices de la celda 4x4
INTEGER :: tipo				! Tipo = Vecindad (viene de Ri)
INTEGER :: J				! Energia de enlace
INTEGER :: i1, i2, i3, j1, j2, j3	! Indices auxiliares
INTEGER :: a, b, c, d, e, f, g, h, k, l, m, n	! Valores auxiliares

! FUNCION
! Inicializacion y energia de enlace
energia = 0.0d0
J = 1

! Entorno para el calculo de las energias
i1 = MODULO(i_cell+1, n_mat)
i2 = MODULO(i_cell-1+n_mat, n_mat)
i3 = MODULO(i_cell+2, n_mat)
j1 = MODULO(j_cell+1, n_mat)
j2 = MODULO(j_cell-1+n_mat, n_mat)
j3 = MODULO(j_cell-2+n_mat, n_mat)

a = mat(i_cell, j_cell)
b = mat(i1, j_cell)
c = mat(i_cell, j2)
d = mat(i1, j2)
e = mat(i_cell, j1)
f = mat(i1, j1)
g = mat(i3,j_cell)
h = mat(i3,j2)
l = mat(i_cell, j3)
k = mat(i1, j3)
n = mat(i2, j_cell)
m = mat(i2, j2)

! Eleccion de la configuracion
SELECT CASE (tipo)
	! Casos con 0 particulas
	CASE (0)	
		a = 0
		b = 0
		c = 0
		d = 0
	! Casos con 1 particula
	CASE (1)	
		a = 1
		b = 0
		c = 0
		d = 0
	CASE (2)	
		a = 0
		b = 1
		c = 0
		d = 0
	CASE (3)	
		a = 0
		b = 0
		c = 1
		d = 0
	CASE (4)	
		a = 0
		b = 0
		c = 0
		d = 1
	! Casos con 2 particulas
	CASE (5)	
		a = 1
		b = 1
		c = 0
		d = 0
	CASE (6)	
		a = 1
		b = 0
		c = 1
		d = 0
	CASE (7)	
		a = 1
		b = 0
		c = 0
		d = 1
	CASE (8)	
		a = 0
		b = 1
		c = 1
		d = 0
	CASE (9)	
		a = 0
		b = 1
		c = 0
		d = 1
	CASE (10)	
		a = 0
		b = 0
		c = 1
		d = 1
	! Casos con 3 particulas
	CASE (11)	
		a = 1
		b = 1
		c = 1
		d = 0
	CASE (12)	
		a = 1
		b = 1
		c = 0
		d = 1
	CASE (13)	
		a = 1
		b = 0
		c = 1
		d = 1
	CASE (14)	
		a = 0
		b = 1
		c = 1
		d = 1
	! Casos con 4 particulas
	CASE (15)	
		a = 1
		b = 1
		c = 1
		d = 1
END SELECT

! Calculo de la energia
energia = energia + KRON(a,b) + KRON(a,c) + KRON(a,e) + KRON(a,n)
energia = energia + KRON(b,f) + KRON(b,g) + KRON(b,d)
energia = energia + KRON(c,m) + KRON(c,d) + KRON(c,l)
energia = energia + KRON(d,h) + KRON(d,k)

! FIN DE LA FUNCION
END FUNCTION ENERGIA


! *********************************************
! 			FUNCTION: E_total
! *********************************************
REAL(8) FUNCTION  E_TOTAL()

! DECLARACION DE VARIABLES LOCALES
implicit none
INTEGER :: i_cell, j_cell 			! Indices para las celdas 2x2
INTEGER :: i1, i2, i3, j1, j2, j3	! Indices para la vecindad
INTEGER :: a, b, c, d, e, f, g, h, k, l, m, n	! Valor de celda vecina

! FUNCION
e_total = 0.0d0
DO i_cell=0,n_mat-1,2
	DO j_cell=0,n_mat-1,2
		
		! Indices para vecinos Von Neumann
		i1 = MODULO(i_cell+1, n_mat)
		i2 = MODULO(i_cell-1+n_mat, n_mat)
		i3 = MODULO(i_cell+2, n_mat)
		j1 = MODULO(j_cell+1, n_mat)
		j2 = MODULO(j_cell-1+n_mat, n_mat)
		j3 = MODULO(j_cell-2+n_mat, n_mat)
		
		! Vecinos Von Neumann
		a = mat(i_cell, j_cell)
		b = mat(i1, j_cell)
		c = mat(i_cell, j2)
		d = mat(i1, j2)
		e = mat(i_cell, j1)
		f = mat(i1, j1)
		g = mat(i3,j_cell)
		h = mat(i3,j2)
		l = mat(i_cell, j3)
		k = mat(i1, j3)
		n = mat(i2, j_cell)
		m = mat(i2, j2)
		
		! Contribucion a la energia total de la configuracion
		e_total = e_total + KRON(a,b) + KRON(a,c)
		e_total = e_total + KRON(b,g) + KRON(b,d)
		e_total = e_total + KRON(c,d) + KRON(c,l)
		e_total = e_total + KRON(d,h) + KRON(d,k)
				
	ENDDO
ENDDO

! FIN DE LA FUNCION
END FUNCTION E_TOTAL


! ====================================================================
! ====================================================================
! FIN DEL MODULO
END MODULE global_var_automata
