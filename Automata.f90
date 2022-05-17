program AUTOMATA
! ====================================================================
! ====================================================================
! Se consideran sistemas con 2 fases:
!	0 -> hueco
!	1 -> materia
!
! Las celdas del sistema evolucionan en grupos de 4 en 4 unidades:
!	|a|b|
!	|c|d|
! Al evolucionar, la celda puede continuar con masa constante o aumentar
! o reducir la masa en una unidad.
!
! Las posibles situaciones de ocupacion se codifican como:
!	codigo = a*1 + b*2 + c*4 + d*8
!
! Se define una regla de evolucion con probabilidades de tipo Arrhenius
! donde la probabilidad de transicion depende de la variacion de energia
! Delta-E = E_final - E_inicial.
!
! La energia de la configuracion se calcula a partir de los vecinos
! de tipo Von Neumann:
!	  * *
!	*|a|b|*
!	*|c|d|*
!	  * *
! ====================================================================
! ====================================================================


! *********************************************
! *********************************************
! 			DECLARACION DE VARIABLES
! *********************************************
! *********************************************
! Se cargan las variables a partir de un modulo
USE global_var_automata
implicit none
INTEGER :: i, j	! Indice para bucle DO


! *********************************************
! *********************************************
! 					MAIN CODE
! *********************************************
! *********************************************
WRITE(*,*) 'Running AUTOMATA'
WRITE(*,*) ' '

! LECTURA DE PARAMETROS DE CONFIGURACION
CALL LECTURA_INPUTS

! DEFINICION DEL DIRECTORIO DE SALVADO
CALL DEFINE_PATH

! INICIALIZACION Y LECTURA DE LA CONFIGURACION INICAL
CALL INICIALIZAR

! EVOLUCION DEL SISTEMA
! Inicializacion de variables
iter = 1
delta_E_rel = 1.0e6

! Calculo de la energia inicial
E_inicial = E_TOTAL()
E_final = E_inicial
E_1 = 1.0d0

! Calculo de la porosidad inicial
porosidad = 0.0d0
DO i=0,n_mat-1,1
	DO j=0,n_mat-1,1
		IF (mat(i,j).EQ.0) porosidad = porosidad + 1
	ENDDO
ENDDO
porosidad = porosidad/(n_mat*n_mat)

! Escritura de valores iniciales
WRITE(21,*) iter, E_inicial
WRITE(26,*) iter, porosidad
WRITE(28,*) iter, E_minima
WRITE(28,*) iter, E_minima	! Repetido para representar la recta
WRITE(29,*) iter, masa
WRITE(30,*) iter, E_inicial/masa

! Bucle principal
DO WHILE (iter.LT.iter_min) !delta_E_rel.GT.crit_end .OR. 
	CALL EVOLUCION			! Evoluciona el sistema
	E_final = E_TOTAL()		! Se recalcula la energia total
	masa = SUM(mat)			! Se recalcula la masa
	
	! Cada iter_prom iteraciones, se promedia la energia
	ener_prom(MODULO(iter,iter_prom)) = E_final
	IF (MODULO(iter,iter_prom).EQ.0) THEN
		E_2 = SUM(ener_prom)/iter_prom
		delta_E_rel = ABS(E_2-E_1)/E_1 * 100
		WRITE(27,*) iter, E_2, delta_E_rel
		E_1 = E_2
	ENDIF
	
	! Actualizacion de la energia y la iteracion
	iter = iter + 1	
	E_inicial = E_final
	
	! Actualizacion de la porosidad
	porosidad = 0.0d0
	DO i=0,n_mat-1,1
		DO j=0,n_mat-1,1
			IF (mat(i,j).EQ.0) porosidad = porosidad + 1
		ENDDO
	ENDDO
	porosidad = porosidad/(n_mat*n_mat)
	
	! Escritura de resultados en fichero
	WRITE(21,*) iter, E_inicial
	WRITE(26,*) iter, porosidad
	WRITE(24,*) iter, mov
	WRITE(29,*) iter, masa
	WRITE(30,*) iter, E_inicial/masa
ENDDO

!~ ! Calculo del promedio temporal
!~ DO i=0,n_mat-1,1
!~ 	DO j=0,n_mat-1,1
!~ 		mat_prom(i,j) = NINT(DBLE(mat_prom(i,j))/iter)
!~ 	ENDDO
!~ 	WRITE(25,*) mat_prom(i,:)
!~ ENDDO

WRITE(*,*) 'La porosidad final es: ', porosidad
WRITE(*,*) 'La energia final es: ', E_inicial
WRITE(*,*) 'La energia minima es: ', E_minima

! CERRADO DE FICHEROS DE ESCRITURA
WRITE(*,*) 'Cerrado de ficheros de escritura'
CLOSE(21)
CLOSE(24)
CLOSE(25)
CLOSE(26)
CLOSE(27)
CLOSE(28)
CLOSE(29)
CLOSE(30)
WRITE(*,*) ' '

! REPRESENTACIONES DEL ESTADO FINAL Y MAGNITUDES
CALL VISUALIZACION_CONF_FINAL
CALL VISUALIZACION_MAGNITUDES
!~ CALL system('gnuplot -e "'// &
!~ 		"m='resultados/promedio-temporal.dat'"// &
!~ 		'" -p gnuplot/volcado.plt')	! Visualizacion con input

WRITE(*,*) '---------------'
WRITE(*,*) 'Done!'



CONTAINS
! *********************************************
! 			SUBROUTINE: Lectura de inputs
! *********************************************
SUBROUTINE LECTURA_INPUTS
USE global_var_automata

! DECLARACION DE VARIABLES LOCALES
implicit none
INTEGER :: i, j	! Indices bucles DO

! SUBRUTINA
! Lectura de inputs
WRITE(*,*) 'Leyendo inputs'
OPEN(UNIT=11, FILE='automata.in', STATUS='OLD')
	READ(11,*) system_filename
	READ(11,*) densidad
	READ(11,*) n_mat
	READ(11,*) T
	READ(11,*) alpha_mas
	READ(11,*) alpha_menos
	READ(11,*) crit_end
	READ(11,*) iter_prom
	READ(11,*) iter_min
	READ(11,*) n_phases
	WRITE(*,*) 'Configuracion inicial: ', system_filename
	WRITE(*,*) 'Tamanio del sistema: ', n_mat
	WRITE(*,*) 'Temperatura: ', T
	WRITE(*,*) 'Coeficiente para aumentar masa: ', alpha_mas
	WRITE(*,*) 'Coeficiente para disminuir masa: ', alpha_menos
	WRITE(*,*) 'Criterio de energia (%) para estado final: ', crit_end
	WRITE(*,*) 'Numero de iteraciones para promedio: ', iter_prom
	WRITE(*,*) 'Numero de iteraciones minimas: ', iter_min
	WRITE(*,*) 'Numero de fases: ', n_phases

	ALLOCATE(kp(0:n_phases-1,0:n_phases-1))
	kp = 0.0d0
	DO i=0,n_phases-1
		DO j=i,n_phases-1
			READ(11,*) kp(i,j)
			kp(j,i) = kp(i,j)
			WRITE(*,*) 'Energia de enlace ',i, '-',j,': ',kp(i,j)
		ENDDO
	ENDDO
CLOSE(11)

! FIN DE LA SUBRUTINA
END SUBROUTINE LECTURA_INPUTS


! *********************************************
! 			SUBROUTINE: Define path
! *********************************************
SUBROUTINE DEFINE_PATH
USE global_var_automata

! DECLARACION DE VARIABLES LOCALES
CHARACTER(len=30) :: T_char, E_char, time_char, &
	alpha_mas_char, alpha_menos_char

! SUBRUTINA
! Pasar variables numéricas a cadena de caracteres
WRITE(T_char,*) T
WRITE(E_char,*) kp(0,1)
WRITE(time_char,*) iter_min
WRITE(alpha_mas_char,'(F9.6)') alpha_mas
WRITE(alpha_menos_char,'(F9.6)') alpha_menos
!WRITE(rho_char,*) densidad

! Crear la carpeta de resultados
CALL getcwd(cwd)
save_dir = 'resultados/' // &
			'E_' // trim(adjustl(E_char(1:9))) // &
			'_alpha_mas_' // trim(adjustl(alpha_mas_char(1:9))) // &
			'_alpha_menos_' // trim(adjustl(alpha_menos_char(1:9))) // &
			'_T_' // trim(adjustl(T_char(1:9))) // &
			'_steps_' // trim(adjustl(time_char(1:20))) 
			
CALL system('mkdir -p ' // save_dir)
CALL SLEEP(1)

! Generacion de ficheros de escritura
WRITE(*,*) 'Inicializacion de ficheros de escritura'
CALL getcwd(cwd)
OPEN(UNIT=21, &
		FILE=trim(cwd) // '/' // trim(save_dir) // &
		'/energia.dat', &
		STATUS='UNKNOWN')
OPEN(UNIT=24, &
		FILE=trim(cwd) // '/' // trim(save_dir) // &
		'/movilidad.dat', &
		STATUS='UNKNOWN')
!~ OPEN(UNIT=25, &
!~ 		FILE=trim(cwd) // '/' // trim(save_dir) // &
!~ 		'/promedio-temporal.dat', &
!~ 		STATUS='UNKNOWN')
OPEN(UNIT=26, &
		FILE=trim(cwd) // '/' // trim(save_dir) // &
		'/porosidad.dat', &
		STATUS='UNKNOWN')
OPEN(UNIT=27, &
		FILE=trim(cwd) // '/' // trim(save_dir) // &
		'/energia-promedio.dat', &
		STATUS='UNKNOWN')
OPEN(UNIT=28, &
		FILE=trim(cwd) // '/' // trim(save_dir) // &
		'/energia-minima.dat', &
		STATUS='UNKNOWN')
OPEN(UNIT=29, &
		FILE=trim(cwd) // '/' // trim(save_dir) // &
		'/masa.dat', &
		STATUS='UNKNOWN')
OPEN(UNIT=30, &
		FILE=trim(cwd) // '/' // trim(save_dir) // &
		'/energia-especifica.dat', &
		STATUS='UNKNOWN')
WRITE(*,*) ' '

! FIN DE LA SUBRUTINA
END SUBROUTINE DEFINE_PATH


! *********************************************
! 	SUBROUTINE: Inicializar
! *********************************************
SUBROUTINE INICIALIZAR
USE global_var_automata

! SUBRUTINA
! Inicializacion de la energia promedio y el estado inicial 
ALLOCATE(ener_prom(0:iter_prom-1))
ener_prom = 0.0d0

ALLOCATE(mat(0:n_mat-1,0:n_mat-1))
ALLOCATE(mat_prom(0:n_mat-1,0:n_mat-1))
mat = 0
mat_prom = 0

IF (system_filename == 'aleatorio') THEN
	! Generacion de una configuracion inicial aleatoria
	WRITE(*,*) 'Configuracion inicial aleatoria, densidad = ', densidad
	CALL CONF_ALEATORIA
	
	! Representacion
	CALL system('gnuplot -e "'// &
		"m='" // TRIM(save_dir) // &
		"/conf-aleatoria-inicial.dat'"// &
		'" -p gnuplot/volcado.plt')	! Visualizacion con input
ELSE
	! Lectura del estado inicial a partir de un fichero de datos
	CALL getcwd(cwd)
	OPEN(UNIT=12, &
			FILE=trim(cwd)//'/sistemas/'//system_filename, &
			STATUS='OLD')
		DO i=0,n_mat-1
			READ(12,*) mat(i,:)
		ENDDO
	CLOSE(12)
	
	! Calculo de la densidad
	densidad = SUM(mat)
	densidad = densidad/(n_mat*n_mat)
	
	! Se almacena una copia de la configuracion inicial en el directorio
	! de resultados
	CALL getcwd(cwd)
	OPEN(UNIT=22, &
			FILE=trim(cwd) // '/' // TRIM(save_dir) // &
			'/config-inicial.dat', &
			STATUS='UNKNOWN')
	DO i=0,n_mat-1
		WRITE(22,*) mat(i,:)
	ENDDO
	CLOSE(22)
	
	!Representacion
	CALL system('gnuplot -e "'// &
		"m='"// TRIM(save_dir) // "/config-inicial.dat'" // &
		'" -p gnuplot/volcado.plt')	! Visualizacion con input
ENDIF
WRITE(*,*) 'Estado inicial cargado'
WRITE(*,*) ' '

! Calculo de la masa del sistema y la energia minima
masa = SUM(mat)
E_minima = 2.0d0 * SQRT(pi) * SQRT(DBLE(masa))

! FIN DE LA SUBRUTINA
END SUBROUTINE INICIALIZAR


! *********************************************
! 	SUBROUTINE: Configuracion aleatoria
! *********************************************
SUBROUTINE CONF_ALEATORIA
USE global_var_automata, ONLY : n_mat, mat, densidad, cwd

! DECLARACION DE VARIABLES LOCALES
implicit none
REAL(8) :: cont
REAL(8) :: i, j	! Nuemros aleatorios
INTEGER :: k	! Indice bucle DO

! SUBRUTINA
! Se genera una configuracion aleatoria de partida
cont = 0
DO WHILE (cont.LT.(n_mat*n_mat*densidad))
	CALL RANDOM_NUMBER(i)	! Numeros aleatorios en [0,1)
	CALL RANDOM_NUMBER(j)
	i = n_mat*i		! Numeros aleatorios entre [0,n_mat)
	j = n_mat*j
	
	! Se rellena la matriz hasta alcanzar la densidad de entrada
	IF (mat(INT(i),INT(j)).EQ.0) THEN
		mat(INT(i),INT(j)) = 1
		cont = cont + 1
	ENDIF
ENDDO

! Se almacena la configuracion inicial en un fichero
CALL getcwd(cwd)
OPEN(UNIT=22, &
		FILE=trim(cwd) // '/' // trim(save_dir) // &
		'/conf-aleatoria-inicial.dat', &
		STATUS='UNKNOWN')
DO k=0,n_mat-1
	WRITE(22,*) mat(k,:)
ENDDO
CLOSE(22)

! FIN DE LA SUBRUTINA
END SUBROUTINE CONF_ALEATORIA

! *********************************************
! 			SUBROUTINE: Visualizacion
! *********************************************
SUBROUTINE VISUALIZACION_CONF_FINAL
USE global_var_automata, ONLY : n_mat, mat, cwd

! DECLARACION DE VARIABLES LOCALES
implicit none
INTEGER :: i! Indice para bucle DO

! SUBRUTINA
! Genera el fichero con la matriz actual
CALL getcwd(cwd)
OPEN(UNIT=23, &
		FILE=trim(cwd) // '/' // trim(save_dir) // &
		'/data-volcado.dat', &
		STATUS='UNKNOWN')
DO i=0,n_mat-1
	WRITE(23,*) mat(i,:)
ENDDO
CLOSE(23)

! Ejecucion de GNUPLOT para representacion
CALL system('gnuplot -e "'// &
		"m='"// TRIM(save_dir) // "/data-volcado.dat'" // &
		'" -p gnuplot/volcado.plt')	! Visualizacion con input

! FIN DE LA SUBRUTINA
END SUBROUTINE VISUALIZACION_CONF_FINAL


! *********************************************
! 		SUBROUTINE: Mag_Visualizacion
! *********************************************
SUBROUTINE VISUALIZACION_MAGNITUDES

! DECLARACION DE VARIABLES LOCALES
implicit none

! SUBRUTINA
! Ejecucion de GNUPLOT para representacion
! Energia y energia promedio
CALL system('gnuplot -e "'// &
		"m='"// TRIM(save_dir) // "/energia.dat'; " // &
		"n='"// TRIM(save_dir) // "/energia-promedio.dat'" // &
		'" -p gnuplot/energia.plt')	

! Movilidad
CALL system('gnuplot -e "'// &
		"m='"// TRIM(save_dir) // "/movilidad.dat'" // &
		'" -p gnuplot/movilidad.plt')

! Masa
CALL system('gnuplot -e "'// &
		"m='"// TRIM(save_dir) // "/masa.dat'" // &
		'" -p gnuplot/masa.plt')

! Energia especifica
CALL system('gnuplot -e "'// &
		"m='"// TRIM(save_dir) // "/energia-especifica.dat'" // &
		'" -p gnuplot/energia-especifica.plt')

! Porosidad
CALL system('gnuplot -e "'// &
		"m='"// TRIM(save_dir) // "/porosidad.dat'" // &
		'" -p gnuplot/porosidad.plt')

! FIN DE LA SUBRUTINA
END SUBROUTINE VISUALIZACION_MAGNITUDES


! *********************************************
! 			SUBROUTINE: Evolucion
! *********************************************
SUBROUTINE EVOLUCION
USE global_var_automata, ONLY : iter, n_mat, mat, kp, T

! DECLARACION DE VARIABLES
implicit none
INTEGER :: i0, j0	! Indices para la distribucion de los grupos de
					! 4x4 celdas que evolucionan
INTEGER :: i, j		! Indices para los bucles DO
INTEGER :: nv		! Numero de ocupacion de la celda 4x4
INTEGER, DIMENSION(:,:),  ALLOCATABLE :: mat_new ! Sistema evolucionado
REAL(8) :: rnd_1, rnd_2		! Numeros aleatorios para crecer materia

! SUBRUTINA
WRITE(*,*) 'Evolucionando, iter = ', iter

! Inicializacion de mat_new y movilidad
ALLOCATE(mat_new(0:n_mat-1,0:n_mat-1))
mat_new = mat
mov = 0.0d0

! Definicion de los grupos de celdas 4x4
!	|A|B|
!	|C|D|
SELECT CASE (MODULO(iter,4))
   CASE (0)	!A
	  i0 = 0
	  j0 = 0
   CASE (1)	!B
	  i0 = 1
	  j0 = 0
   CASE (2)	!D
	  i0 = 1
	  j0 = 1
   CASE (3)	!C
	  i0 = 0
	  j0 = 1
END SELECT

! Bucle principal para recorrer todas las celdad 4x4
DO i=i0,n_mat-1,2
	DO j=j0,n_mat-1,2
		! Numero de ocupacion de la celda 4x4
		nv = mat(i,j) + &
			mat(MODULO(i+1,n_mat),j) + &
			mat(i,MODULO(j-1+n_mat,n_mat)) + &
			mat(MODULO(i+1,n_mat),MODULO(j-1+n_mat,n_mat))
			
		! Evolucion segun el numero de ocupacion
		SELECT CASE (nv)
			CASE (0)	
				CALL R0(i, j, mat_new)
			CASE (1)
				CALL R1(i, j, mat_new)
			CASE (2)
				CALL R2(i, j, mat_new)
			CASE (3)
				CALL R3(i, j, mat_new)
			CASE (4)	
				CALL R4(i, j, mat_new)
		END SELECT
	ENDDO
ENDDO

! Calculo de la movilidad y el promedio temporal para la porosidad
DO i=0,n_mat-1,1
	DO j=0,n_mat-1,1
		mov = mov + MODULO(mat(i,j)+mat_new(i,j), 2)
		IF (mat_new(i,j).EQ.1) mat_prom(i,j) = mat_prom(i,j) + 1
	ENDDO
ENDDO
mov = mov / (n_mat*n_mat)

! Actualizacion de la matriz
mat = mat_new

!~ ! FIN DE LA SUBRUTINA
END SUBROUTINE EVOLUCION


! *********************************************
! 				SUBROUTINE: R0
! *********************************************
SUBROUTINE R0(i_cell, j_cell, mat_fin)
USE global_var_automata, ONLY : n_mat, mat, alpha_mas

! DECLARACION DE VARIABLES
implicit none
INTEGER :: i_cell, j_cell	! Indices de la celda 4x4
INTEGER, DIMENSION(:,:),  ALLOCATABLE :: mat_fin ! Sistema final
INTEGER :: vecindad		! Vecindad
REAL(8) :: E_ini_cell	! Energia de la configuracion inicial
REAL(8) :: E_fin_cell	! Energia de la configuracion final
REAL(8), DIMENSION(0:4) :: D_ener_tran_k ! Energia de la transicion k
REAL(8), DIMENSION(0:4) :: prob_tran_k ! Probabil. de la transicion k
REAL(8) :: pt	! Probabilidad total para normalizacion
REAL(8), DIMENSION(0:4) :: q_k	! Probabilidades normalizadas
INTEGER :: a, b, c, d	! Fase en cada punto de la celda 4x4
INTEGER :: k	! Indice para bucle DO 
REAL(8) :: aux	! Numero auxiliar para la variable aleatoria
INTEGER, DIMENSION(1:4) :: estados	! Estados posibles al aumentar masa

! SUBRUTINA
! Vecindad inicial y estados posibles si aumenta la masa
vecindad = 0
estados = (/1,2,3,4/)

! Calculo de la energia inicial 
E_ini_cell = ENERGIA(i_cell, j_cell, vecindad)

! Probabilidad para cada transicion con masa constante
D_ener_tran_k = 0
prob_tran_k = 0
pt = 0

E_fin_cell = ENERGIA(i_cell, j_cell, vecindad)
D_ener_tran_k(0) = E_fin_cell - E_ini_cell
prob_tran_k(0) = exp(-D_ener_tran_k(0)/T)
pt = pt + prob_tran_k(0)

! Probabilidades para cada transicion con aumento de masa
! Son accesibles los estados con codigo {1,2,4,8}
DO k=1,4
	E_fin_cell = ENERGIA(i_cell, j_cell, estados(k))
	D_ener_tran_k(k) = E_fin_cell - E_ini_cell
	prob_tran_k(k) = alpha_mas * exp(-D_ener_tran_k(k)/T)
	pt = pt + prob_tran_k(k)
ENDDO

! Normalizacion de probabilidades
DO k=0,4
	prob_tran_k(k) = prob_tran_k(k)  / pt
ENDDO

q_k(0) = prob_tran_k(0)
DO k=1,4
	q_k(k) = q_k(k-1) + prob_tran_k(k)
ENDDO

! Eleccion de la transicion de acuerdo a su probabilidad
CALL RANDOM_NUMBER(aux)
IF (aux.LE.q_k(0)) THEN
	a = 0
	b = 0
	c = 0
	d = 0
ELSEIF (aux.GT.q_k(0) .AND. aux.LE.q_k(1)) THEN
	a = 1
	b = 0
	c = 0
	d = 0
ELSEIF (aux.GT.q_k(1) .AND. aux.LE.q_k(2)) THEN
	a = 0
	b = 1
	c = 0
	d = 0
ELSEIF (aux.GT.q_k(2) .AND. aux.LE.q_k(3)) THEN
	a = 0
	b = 0
	c = 1
	d = 0
ELSE
	a = 0
	b = 0
	c = 0
	d = 1
ENDIF

! Actualizacion de la celda tras evolucionar
mat_fin(i_cell,j_cell) = a
mat_fin(MODULO(i_cell+1,n_mat),j_cell) = b
mat_fin(i_cell,MODULO(j_cell-1+n_mat,n_mat)) = c
mat_fin(MODULO(i_cell+1,n_mat),MODULO(j_cell-1+n_mat,n_mat)) = d

! FIN DE LA SUBRUTINA
END SUBROUTINE R0


! *********************************************
! 				SUBROUTINE: R1
! *********************************************
SUBROUTINE R1(i_cell, j_cell, mat_fin)
USE global_var_automata, ONLY : n_mat, mat, ENERGIA, &
									alpha_mas, alpha_menos

! DECLARACION DE VARIABLES
implicit none
INTEGER :: i_cell, j_cell	! Indices de la celda 4x4
INTEGER, DIMENSION(:,:),  ALLOCATABLE :: mat_fin ! Sistema final
INTEGER :: codigo		! Codigo para la vecindad
INTEGER :: vecindad		! Vecindad
REAL(8) :: E_ini_cell	! Energia de la configuracion inicial
REAL(8) :: E_fin_cell	! Energia de la configuracion final
REAL(8), DIMENSION(1:8) :: D_ener_tran_k ! Energia de la transicion k
REAL(8), DIMENSION(1:8) :: prob_tran_k ! Probabil. de la transicion k
REAL(8) :: pt	! Probabilidad total para normalizacion
REAL(8), DIMENSION(1:8) :: q_k	! Probabilidades normalizadas
INTEGER :: a, b, c, d	! Fase en cada punto de la celda 4x4
INTEGER :: k	! Indice para bucle DO 
REAL(8) :: aux	! Numero auxiliar para la variable aleatoria
INTEGER, DIMENSION(1:3) :: estados_mas! Estados posibles al crecer masa
INTEGER, DIMENSION(1:1) :: estados_menos! Estados posibles al bajar masa

! SUBRUTINA
! Codigo para vecindades: {1,2,4,8}
codigo = mat(i_cell, j_cell) + &
	mat(MODULO(i_cell+1,n_mat), j_cell) * 2 + &
	mat(i_cell, MODULO(j_cell-1+n_mat,n_mat)) * 4 + &
	mat(MODULO(i_cell+1,n_mat), MODULO(j_cell-1+n_mat,n_mat)) * 8

! Definicion de vecindad y estados con variacion de masa
SELECT CASE (codigo)
	CASE (1)
		vecindad = 1
		estados_mas = (/5,6,7/)
	CASE (2)
		vecindad = 2
		estados_mas = (/5,8,9/)
	CASE (4)
		vecindad = 3
		estados_mas = (/6,8,10/)
	CASE (8)
		vecindad = 4
		estados_mas = (/7,9,10/)
	CASE DEFAULT
		WRITE(*,*) 'Error en SUBRUTINA R1: codigo no valido'
END SELECT
estados_menos = (/0/)

! Calculo de la energia inicial 
E_ini_cell = ENERGIA(i_cell, j_cell, vecindad)

! Probabilidad para cada transicion con masa constante
D_ener_tran_k = 0
prob_tran_k = 0
pt = 0

DO k=1,4
	E_fin_cell = ENERGIA(i_cell, j_cell, k)
	D_ener_tran_k(k) = E_fin_cell - E_ini_cell
	prob_tran_k(k) = exp(-D_ener_tran_k(k)/T)
	pt = pt + prob_tran_k(k)
ENDDO

! Probabilidades para cada transicion con aumento de masa
! Son accesibles los estados con codigo {3,5,6,9,10,12}
DO k=1,3
	E_fin_cell = ENERGIA(i_cell, j_cell, estados_mas(k))
	D_ener_tran_k(4+k) = E_fin_cell - E_ini_cell
	prob_tran_k(4+k) = alpha_mas * exp(-D_ener_tran_k(4+k)/T)
	pt = pt + prob_tran_k(4+k)
ENDDO

! Probabilidades para cada transicion con disminucion de masa
! Son accesibles los estados con codigo {0}
E_fin_cell = ENERGIA(i_cell, j_cell, estados_menos(1))
D_ener_tran_k(8) = E_fin_cell - E_ini_cell
prob_tran_k(8) = alpha_menos * exp(-D_ener_tran_k(8)/T)
pt = pt + prob_tran_k(8)

! Normalizacion de probabilidades
DO k=1,8
	prob_tran_k(k) = prob_tran_k(k)  / pt
ENDDO

q_k(1) = prob_tran_k(1)
DO k=2,8
	q_k(k) = q_k(k-1) + prob_tran_k(k)
ENDDO

! Eleccion de la transicion de acuerdo a su probabilidad
CALL RANDOM_NUMBER(aux)
IF (aux.LE.q_k(1)) THEN
	a = 1
	b = 0
	c = 0
	d = 0
ELSEIF (aux.GT.q_k(1) .AND. aux.LE.q_k(2)) THEN
	a = 0
	b = 1
	c = 0
	d = 0
ELSEIF (aux.GT.q_k(2) .AND. aux.LE.q_k(3)) THEN
	a = 0
	b = 0
	c = 1
	d = 0
ELSEIF (aux.GT.q_k(3) .AND. aux.LE.q_k(4)) THEN
	a = 0
	b = 0
	c = 0
	d = 1
ELSEIF (aux.GT.q_k(4) .AND. aux.LE.q_k(5)) THEN
	vecindad = estados_mas(1)
ELSEIF (aux.GT.q_k(5) .AND. aux.LE.q_k(6)) THEN
	vecindad = estados_mas(2)
ELSEIF (aux.GT.q_k(6) .AND. aux.LE.q_k(7)) THEN
	vecindad = estados_mas(3)
ELSE
	vecindad = estados_menos(1)
ENDIF

! Configuraciones finales si se aumenta o disminuye la masa
SELECT CASE (vecindad)
	CASE(5)
		a = 1
		b = 1
		c = 0
		d = 0
	CASE(6)
		a = 1
		b = 0
		c = 1
		d = 0
	CASE(7)
		a = 1
		b = 0
		c = 0
		d = 1
	CASE(8)
		a = 0
		b = 1
		c = 1
		d = 0
	CASE(9)
		a = 0
		b = 1
		c = 0
		d = 1
	CASE(10)
		a = 0
		b = 0
		c = 1
		d = 1
	CASE(0)
		a = 0
		b = 0
		c = 0
		d = 0
END SELECT

! Actualizacion de la celda  4x4 tras evolucionar
mat_fin(i_cell, j_cell) = a
mat_fin(MODULO(i_cell+1,n_mat), j_cell) = b
mat_fin(i_cell, MODULO(j_cell-1+n_mat,n_mat)) = c
mat_fin(MODULO(i_cell+1,n_mat), MODULO(j_cell-1+n_mat,n_mat)) = d

! FIN DE LA SUBRUTINA
END SUBROUTINE R1


! *********************************************
! 				SUBROUTINE: R2
! *********************************************
SUBROUTINE R2(i_cell, j_cell, mat_fin)
USE global_var_automata, ONLY : n_mat, mat, ENERGIA, &
									alpha_mas, alpha_menos

! DECLARACION DE VARIABLES
implicit none
INTEGER :: i_cell, j_cell	! Indices de la celda 4x4
INTEGER, DIMENSION(:,:),  ALLOCATABLE :: mat_fin ! Sistema final
INTEGER :: codigo		! Codigo para la vecindad
INTEGER :: vecindad		! Vecindad
REAL(8) :: E_ini_cell	! Energia de la configuracion inicial
REAL(8) :: E_fin_cell	! Energia de la configuracion final
REAL(8), DIMENSION(5:14) :: D_ener_tran_k	! Energia de la transicion k
REAL(8), DIMENSION(5:14) :: prob_tran_k ! Probabil. de la transicion k
REAL(8) :: pt	! Probabilidad total para normalizacion
REAL(8), DIMENSION(5:14) :: q_k	! Probabilidades normalizadas
INTEGER :: a, b, c, d	! Fase en cada punto de la celda 4x4
INTEGER :: k	! Indice para bucle DO 
REAL(8) :: aux	! Numero auxiliar para la variable aleatoria
INTEGER, DIMENSION(1:2) :: estados_mas! Estados posibles al crecer masa
INTEGER, DIMENSION(1:2) :: estados_menos! Estados posibles al bajar masa

! SUBRUTINA
! Codigo para vecindades: {3,5,6,9,10,12}
codigo = mat(i_cell, j_cell) + &
	mat(MODULO(i_cell+1,n_mat), j_cell) * 2 + &
	mat(i_cell, MODULO(j_cell-1+n_mat,n_mat)) * 4 + &
	mat(MODULO(i_cell+1,n_mat), MODULO(j_cell-1+n_mat,n_mat)) * 8

! Definicion de vecindad
SELECT CASE(codigo)
	CASE(3)
		vecindad = 5
		estados_mas = (/11,12/)
		estados_menos = (/1,2/)
	CASE(5)
		vecindad = 6
		estados_mas = (/11,13/)
		estados_menos = (/1,3/)
	CASE(6)
		vecindad = 8
		estados_mas = (/11,14/)
		estados_menos = (/2,3/)
	CASE(9)
		vecindad = 7
		estados_mas = (/12,13/)
		estados_menos = (/1,4/)
	CASE(10)
		vecindad = 9
		estados_mas = (/12,14/)
		estados_menos = (/2,4/)
	CASE(12)
		vecindad = 10
		estados_mas = (/13,14/)
		estados_menos = (/3,4/)
	CASE DEFAULT
		WRITE(*,*) 'Error en SUBRUTINA R2: codigo no valido'
END SELECT

! Calculo de la energia inicial 
E_ini_cell = ENERGIA(i_cell, j_cell, vecindad)

! Probabilidad para cada transicion con masa constante
D_ener_tran_k = 0
prob_tran_k = 0
pt = 0

DO k=5,10
	E_fin_cell = ENERGIA(i_cell, j_cell, k)
	D_ener_tran_k(k) = E_fin_cell - E_ini_cell
	prob_tran_k(k) = exp(-D_ener_tran_k(k)/T)
	pt = pt + prob_tran_k(k)
ENDDO

! Probabilidades para cada transicion con aumento de masa
! Son accesibles los estados con codigo {7,11,13,14}
DO k=1,2
	E_fin_cell = ENERGIA(i_cell, j_cell, estados_mas(k))
	D_ener_tran_k(10+k) = E_fin_cell - E_ini_cell
	prob_tran_k(10+k) = alpha_mas * exp(-D_ener_tran_k(10+k)/T)
	pt = pt + prob_tran_k(10+k)
ENDDO

! Probabilidades para cada transicion con reduccion de masa
! Son accesibles los estados con codigo {1,2,4,8}
DO k=1,2
	E_fin_cell = ENERGIA(i_cell, j_cell, estados_menos(k))
	D_ener_tran_k(12+k) = E_fin_cell - E_ini_cell
	prob_tran_k(12+k) = alpha_menos * exp(-D_ener_tran_k(12+k)/T)
	pt = pt + prob_tran_k(12+k)
ENDDO

! Normalizacion de probabilidades
DO k=5,14
	prob_tran_k(k) = prob_tran_k(k) / pt
ENDDO

q_k(5) = prob_tran_k(5)
DO k=6,14
	q_k(k) = q_k(k-1) + prob_tran_k(k)
ENDDO

! Eleccion de la transicion de acuerdo a su probabilidad
CALL RANDOM_NUMBER(aux)
IF (aux.LE.q_k(5)) THEN
	a = 1
	b = 1
	c = 0
	d = 0
ELSEIF (aux.GT.q_k(5) .AND. aux.LE.q_k(6)) THEN
	a = 1
	b = 0
	c = 1
	d = 0
ELSEIF (aux.GT.q_k(6) .AND. aux.LE.q_k(7)) THEN
	a = 1
	b = 0
	c = 0
	d = 1
ELSEIF (aux.GT.q_k(7) .AND. aux.LE.q_k(8)) THEN
	a = 0
	b = 1
	c = 1
	d = 0
ELSEIF (aux.GT.q_k(8) .AND. aux.LE.q_k(9)) THEN
	a = 0
	b = 1
	c = 0
	d = 1
ELSEIF (aux.GT.q_k(9) .AND. aux.LE.q_k(10)) THEN
	a = 0
	b = 0
	c = 1
	d = 1
ELSEIF (aux.GT.q_k(10) .AND. aux.LE.q_k(11)) THEN
	vecindad = estados_mas(1)
ELSEIF (aux.GT.q_k(11) .AND. aux.LE.q_k(12)) THEN
	vecindad = estados_mas(2)
ELSEIF (aux.GT.q_k(12) .AND. aux.LE.q_k(13)) THEN
	vecindad = estados_menos(1)
ELSE
	vecindad = estados_menos(2)
ENDIF

! Configuraciones finales si se aumenta o disminuye la masa
SELECT CASE (vecindad)
	CASE(11)
		a = 1
		b = 1
		c = 1
		d = 0
	CASE(12)
		a = 1
		b = 1
		c = 0
		d = 1
	CASE(13)
		a = 1
		b = 0
		c = 1
		d = 1
	CASE(14)
		a = 0
		b = 1
		c = 1
		d = 1
	CASE(1)
		a = 1
		b = 0
		c = 0
		d = 0
	CASE(2)
		a = 0
		b = 1
		c = 0
		d = 0
	CASE(3)
		a = 0
		b = 0
		c = 1
		d = 0
	CASE(4)
		a = 0
		b = 0
		c = 0
		d = 1
END SELECT

! Actualizacion de la celda  4x4 tras evolucionar
mat_fin(i_cell, j_cell) = a
mat_fin(MODULO(i_cell+1,n_mat), j_cell) = b
mat_fin(i_cell, MODULO(j_cell-1+n_mat,n_mat)) = c
mat_fin(MODULO(i_cell+1,n_mat), MODULO(j_cell-1+n_mat,n_mat)) = d

! FIN DE LA SUBRUTINA
END SUBROUTINE R2


! *********************************************
! 				SUBROUTINE: R3
! *********************************************
SUBROUTINE R3(i_cell, j_cell, mat_fin)
USE global_var_automata, ONLY : n_mat, mat, ENERGIA, &
									alpha_mas, alpha_menos

! DECLARACION DE VARIABLES
implicit none
INTEGER :: i_cell, j_cell	! Indices de la celda 4x4
INTEGER, DIMENSION(:,:),  ALLOCATABLE :: mat_fin ! Sistema final
INTEGER :: codigo		! Codigo para la vecindad
INTEGER :: vecindad		! Vecindad
REAL(8) :: E_ini_cell	! Energia de la configuracion inicial
REAL(8) :: E_fin_cell	! Energia de la configuracion final
REAL(8), DIMENSION(11:18) :: D_ener_tran_k	! Energia de la transicion k
REAL(8), DIMENSION(11:18) :: prob_tran_k ! Probabil. de la transicion k
REAL(8) :: pt	! Probabilidad total para normalizacion
REAL(8), DIMENSION(11:18) :: q_k	! Probabilidades normalizadas
INTEGER :: a, b, c, d	! Fase en cada punto de la celda 4x4
INTEGER :: k	! Indice para bucle DO 
REAL(8) :: aux	! Numero auxiliar para la variable aleatoria
INTEGER, DIMENSION(1:1) :: estados_mas! Estados posibles al crecer masa
INTEGER, DIMENSION(1:3) :: estados_menos! Estados posibles al bajar masa
						
! SUBRUTINA
! Codigo para vecindades: {7,11,13,14}
codigo = mat(i_cell, j_cell) + &
	mat(MODULO(i_cell+1,n_mat), j_cell) * 2 + &
	mat(i_cell, MODULO(j_cell-1+n_mat,n_mat)) * 4 + &
	mat(MODULO(i_cell+1,n_mat), MODULO(j_cell-1+n_mat,n_mat)) * 8

! Definicion de vecindad
estados_mas = 15
SELECT CASE (codigo)
	CASE (7)
		vecindad = 11
		estados_menos = (/8,6,5/)
	CASE (11)
		vecindad = 12
		estados_menos = (/9,7,5/)
	CASE (13)
		vecindad = 13
		estados_menos = (/10,7,6/)
	CASE (14)
		vecindad = 14
		estados_menos = (/10,9,8/)
	CASE DEFAULT
		WRITE(*,*) 'Error en SUBRUTINA R3: codigo no valido'
END SELECT

! Calculo de la energia inicial 
E_ini_cell = ENERGIA(i_cell, j_cell, vecindad)

! Probabilidad para cada transicion con masa constante
D_ener_tran_k = 0
prob_tran_k = 0
pt = 0

DO k=11,14
	E_fin_cell = ENERGIA(i_cell, j_cell, k)
	D_ener_tran_k(k) = E_fin_cell - E_ini_cell
	prob_tran_k(k) = exp(-D_ener_tran_k(k)/T)
	pt = pt + prob_tran_k(k)
ENDDO

! Probabilidades para cada transicion con aumento de masa
! Son accesibles los estados con codigo {15}
E_fin_cell = ENERGIA(i_cell, j_cell, estados_mas(1))
D_ener_tran_k(15) = E_fin_cell - E_ini_cell
prob_tran_k(15) = alpha_mas * exp(-D_ener_tran_k(15)/T)
pt = pt + prob_tran_k(15)

! Probabilidades para cada transicion con reduccion de masa
! Son accesibles los estados con codigo {3,5,6,9,10,12}
DO k=1,3
	E_fin_cell = ENERGIA(i_cell, j_cell, estados_menos(k))
	D_ener_tran_k(15+k) = E_fin_cell - E_ini_cell
	prob_tran_k(15+k) = alpha_menos * exp(-D_ener_tran_k(15+k)/T)
	pt = pt + prob_tran_k(15+k)
ENDDO

! Normalizacion de probabilidades
DO k=11,18
	prob_tran_k(k) = prob_tran_k(k) / pt
ENDDO

q_k(11) = prob_tran_k(11)
DO k=12,18
	q_k(k) = q_k(k-1) + prob_tran_k(k)
ENDDO

! Eleccion de la transicion de acuerdo a su probabilidad
CALL RANDOM_NUMBER(aux)
IF (aux.LE.q_k(11)) THEN
	a = 1
	b = 1
	c = 1
	d = 0
ELSEIF (aux.GT.q_k(11) .AND. aux.LE.q_k(12)) THEN
	a = 1
	b = 1
	c = 0
	d = 1
ELSEIF (aux.GT.q_k(12) .AND. aux.LE.q_k(13)) THEN
	a = 1
	b = 0
	c = 1
	d = 1
ELSEIF (aux.GT.q_k(13) .AND. aux.LE.q_k(14)) THEN
	a = 0
	b = 1
	c = 1
	d = 1
ELSEIF (aux.GT.q_k(14) .AND. aux.LE.q_k(15)) THEN
	vecindad = estados_mas(1)
ELSEIF (aux.GT.q_k(15) .AND. aux.LE.q_k(16)) THEN
	vecindad = estados_menos(1)
ELSEIF (aux.GT.q_k(16) .AND. aux.LE.q_k(17)) THEN
	vecindad = estados_menos(2)
ELSE
	vecindad = estados_menos(3)
ENDIF

! Configuraciones finales si se aumenta o disminuye la masa
SELECT CASE (vecindad)
	CASE(15)
		a = 1
		b = 1
		c = 1
		d = 1
	CASE(5)
		a = 1
		b = 1
		c = 0
		d = 0
	CASE(6)
		a = 1
		b = 0
		c = 1
		d = 0
	CASE(7)
		a = 1
		b = 0
		c = 0
		d = 1
	CASE(8)
		a = 0
		b = 1
		c = 1
		d = 0
	CASE(9)
		a = 0
		b = 1
		c = 0
		d = 1
	CASE(10)
		a = 0
		b = 0
		c = 1
		d = 1
END SELECT

! Actualizacion de la celda  4x4 tras evolucionar
mat_fin(i_cell, j_cell) = a
mat_fin(MODULO(i_cell+1,n_mat), j_cell) = b
mat_fin(i_cell, MODULO(j_cell-1+n_mat,n_mat)) = c
mat_fin(MODULO(i_cell+1,n_mat), MODULO(j_cell-1+n_mat,n_mat)) = d

! FIN DE LA SUBRUTINA
END SUBROUTINE R3


! *********************************************
! 				SUBROUTINE: R4
! *********************************************
SUBROUTINE R4(i_cell, j_cell, mat_fin)
USE global_var_automata, ONLY : n_mat, mat, alpha_menos

! DECLARACION DE VARIABLES
implicit none
INTEGER :: i_cell, j_cell	! Indices de la celda 4x4
INTEGER, DIMENSION(:,:),  ALLOCATABLE :: mat_fin ! Sistema final
INTEGER :: codigo		! Codigo para la vecindad
INTEGER :: vecindad		! Vecindad
REAL(8) :: E_ini_cell	! Energia de la configuracion inicial
REAL(8) :: E_fin_cell	! Energia de la configuracion final
REAL(8), DIMENSION(15:19) :: D_ener_tran_k	! Energia de la transicion k
REAL(8), DIMENSION(15:19) :: prob_tran_k ! Probabil. de la transicion k
REAL(8) :: pt	! Probabilidad total para normalizacion
REAL(8), DIMENSION(15:19) :: q_k	! Probabilidades normalizadas
INTEGER :: a, b, c, d	! Fase en cada punto de la celda 4x4
INTEGER :: k	! Indice para bucle DO 
REAL(8) :: aux	! Numero auxiliar para la variable aleatoria
INTEGER, DIMENSION(1:4) :: estados_menos! Estados posibles al bajar masa

! SUBRUTINA
! Vecindad inicial y estados posibles si aumenta la masa
vecindad = 15
estados_menos = (/11,12,13,14/)

! Calculo de la energia inicial 
E_ini_cell = ENERGIA(i_cell, j_cell, vecindad)

! Probabilidad para cada transicion con masa constante
D_ener_tran_k = 0
prob_tran_k = 0
pt = 0

E_fin_cell = ENERGIA(i_cell, j_cell, vecindad)
D_ener_tran_k(15) = E_fin_cell - E_ini_cell
prob_tran_k(15) = exp(-D_ener_tran_k(15)/T)
pt = pt + prob_tran_k(15)

! Probabilidades para cada transicion con reduccion de masa
! Son accesibles los estados con codigo {7,11,13,14}
DO k=1,4
	E_fin_cell = ENERGIA(i_cell, j_cell, estados_menos(k))
	D_ener_tran_k(15+k) = E_fin_cell - E_ini_cell
	prob_tran_k(15+k) = alpha_menos * exp(-D_ener_tran_k(15+k)/T)
	pt = pt + prob_tran_k(15+k)
ENDDO

! Normalizacion de probabilidades
DO k=15,19
	prob_tran_k(k) = prob_tran_k(k)  / pt
ENDDO

q_k(15) = prob_tran_k(15)
DO k=16,19
	q_k(k) = q_k(k-1) + prob_tran_k(k)
ENDDO

! Eleccion de la transicion de acuerdo a su probabilidad
CALL RANDOM_NUMBER(aux)
IF (aux.LE.q_k(15)) THEN
	a = 1
	b = 1
	c = 1
	d = 1
ELSEIF (aux.GT.q_k(15) .AND. aux.LE.q_k(16)) THEN
	a = 1
	b = 1
	c = 1
	d = 0
ELSEIF (aux.GT.q_k(16) .AND. aux.LE.q_k(17)) THEN
	a = 1
	b = 1
	c = 0
	d = 1
ELSEIF (aux.GT.q_k(17) .AND. aux.LE.q_k(18)) THEN
	a = 1
	b = 0
	c = 1
	d = 1
ELSE
	a = 0
	b = 1
	c = 1
	d = 1
ENDIF

! Actualizacion de la celda tras evolucionar
mat_fin(i_cell,j_cell) = a
mat_fin(MODULO(i_cell+1,n_mat),j_cell) = b
mat_fin(i_cell,MODULO(j_cell-1+n_mat,n_mat)) = c
mat_fin(MODULO(i_cell+1,n_mat),MODULO(j_cell-1+n_mat,n_mat)) = d

! FIN DE LA SUBRUTINA
END SUBROUTINE R4


! ====================================================================
! ====================================================================
! FIN DEL PROGRAMA PRINCIPAL
end program AUTOMATA
