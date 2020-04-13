! 18 October 2012 - Lorenzo Lodi
! This simple program, given the spectroscopic term symbols 
!  of two atoms (in L.S or Russel-Sanders coupling) 
!  prints the possible symmetries of the molecular terms of the resulting
!  diatomic molecule which correlate adiabatically with them.
! This is an application of the so-called Wigner-Witmer or correlations
! rules; ! see R.N. Zare, Angular Momentum, p.66-69.
! An atomic term symbol specified by multiplicity (=2S+1), value of L 
! (S=0, P=1, D=2, etc ) and the inversion parity, even or odd.
! Term symbols for all atoms can be obtained online from the 
! NIST Atomic Spectra Database http://www.nist.gov/pml/data/asd.cfm
! (choose then "Levels" and remember that, e.g., Fe I indicate 
!  neutral iron, Fe II indicates Fe^+, etc.
! NOTE: The convention for atomic term symbols is that odd-parity
!       terms have an "o" superscript, while even-parity terms
!       have no superscript at all.
!
! ** EXAMPLE **
! Consider the OH radical; let us see what molecular terms correlated
! with the ground states of the two atoms.
! Hydrogen has a 2 S ground state (even parity)
! Oxygen   has a 3 P ground state (even parity)
! Entering in the program "2 S even" and "3 P even" 
! the program says we get two possible multiplicities ( 2 or 4) and
! one Σ- and one Π state for each multiplicity.
! So overall we have four molecular terms: 2Σ-, 4Σ-, 2Π and 4Π.
!
! The "total number of microstates" number is the total geneneracy.
! When spin-orbit is included some states split. 
! Adding an external electric field splits all degeneracies.
! Transition metals leads to lots of microstates. 
! For example, for the iron diatomic Fe_2 no less than 75 molecular
! curves correlated with the atomic ground state, and the total
! degeneracy is 625! 
!
! To do list:
! 2) g/u labels for homonuclear molecules
! 3) Rules for the spin-orbit split components
! Hydrogen has a ground state with term symbol 1 S
!
program WigmerWitmer
implicit none
integer, parameter :: dp=kind(1.d0)

character(len=1) text_L1, text_L2 ! text representing S, P, D, F...
character(len=1) text_p1, text_p2 ! parity, e=even, o=odd
character(len=1) text_homonuclear ! yes if homonuclear, no if not
integer :: L1, L2 ! atomic angular momenta
integer :: p1, p2 ! atomic parities (+1 or -1)
real(dp) :: S1, S2 ! atomic spin
real(dp) :: Smin, Smax
integer :: M1, M2              ! atomic multiplicities
integer :: g1, g2              ! atomic degeneracies
integer :: nTotal, nLambda, lambda, nLambdaPlus, nLambdaMinus
integer :: i, l, m, kk
integer :: iHomoNuclear

write(*,'(A)') 'October 2012, Lorenzo Lodi'
write(*,'(A)') 'Wigner-Witmer rules for diatomics'
write(*,'(A)') 'TEST: 4Σ+g, 7Πu, Δ'
write(*,'(A)') 'TEST: 4Sigma+g, 7Pi_u, Δ'
write(*,'(A)') 'Enter atomic term symbols: multiplicity & term label (eg, 2 S odd )'

write(*,'(A)',advance='no') 'Atom 1 = '
read(*,*) M1, text_L1, text_p1

write(*,'(A)',advance='no') 'Atom 2 = '
read(*,*) M2, text_L2, text_p2

if(M1 < 0) stop 'FATAL: negative M1'
if(M2 < 0) stop 'FATAL: negative M1'

S1 = real(M1-1,dp)/2._dp
S2 = real(M2-1,dp)/2._dp

L1 = -1
if( text_L1 == 'S' .or. text_L1 == 's') L1 = 0
if( text_L1 == 'P' .or. text_L1 == 'p') L1 = 1
if( text_L1 == 'D' .or. text_L1 == 'd') L1 = 2
if( text_L1 == 'F' .or. text_L1 == 'f') L1 = 3
if( text_L1 == 'G' .or. text_L1 == 'g') L1 = 4
if( text_L1 == 'H' .or. text_L1 == 'h') L1 = 5
if( text_L1 == 'I' .or. text_L1 == 'i') L1 = 6
if( text_L1 == 'K' .or. text_L1 == 'k') L1 = 7
if( text_L1 == 'L' .or. text_L1 == 'l') L1 = 8
if( L1 < 0) stop 'FATAL: I do not like L1 '

L2 = -1
if( text_L2 == 'S' .or. text_L2 == 's') L2 = 0
if( text_L2 == 'P' .or. text_L2 == 'p') L2 = 1
if( text_L2 == 'D' .or. text_L2 == 'd') L2 = 2
if( text_L2 == 'F' .or. text_L2 == 'f') L2 = 3
if( text_L2 == 'G' .or. text_L2 == 'g') L2 = 4
if( text_L2 == 'H' .or. text_L2 == 'h') L2 = 5
if( text_L2 == 'I' .or. text_L2 == 'i') L2 = 6
if( text_L2 == 'K' .or. text_L2 == 'k') L2 = 7
if( text_L2 == 'L' .or. text_L2 == 'l') L2 = 8
if( L2 < 0._dp) stop 'FATAL: I do not like L2 '


p1 = -2
if( text_p1 == 'E' .or. text_p1 == 'e') p1 = +1
if( text_p1 == 'O' .or. text_p1 == 'o') p1 = -1
if( abs(p1) /= 1 ) stop 'FATAL: I do not like p1'


p2 = -2
if( text_p2 == 'E' .or. text_p2 == 'e') p2 = +1
if( text_p2 == 'O' .or. text_p2 == 'o') p2 = -1
if( abs(p2) /= 1 ) stop 'FATAL: I do not like p2'


write(*,'(A)', advance='no') 'Atoms are the same? (yes/no) '
read(*,*) text_homonuclear


iHomoNuclear=-2
if( text_Homonuclear == 'y' .or. text_Homonuclear == 'Y') iHomoNuclear = +1
if( text_Homonuclear == 'n' .or. text_Homonuclear == 'N') iHomoNuclear =  0
if( iHomoNuclear ==-2) stop 'FATAL: I do not like your reply '

g1 = (2*L1+1)*M1
g2 = (2*L2+1)*M2

write(*,*)
write(*,'(A, I5, 4x, A, F5.1,4x,a,I6,a,I3)') 'Atom 1 ::  L = ', L1, ' S = ', S1, 'degeneracy g1 = ', g1, ' parity =', p1
write(*,'(A, I5, 4x, A, F5.1,4x,a,I6,a,I3)') 'Atom 2 ::  L = ', L2, ' S = ', S2, 'degeneracy g2 = ', g2, ' parity =', p2
write(*,*)
write(*,'(A,I6)') 'Total number of microstates g1*g2 = ', g1*g2
write(*,*)
Smin = min(S1, S2)
Smax = max(S1, S2)

nTotal = min(M1,M2)*( 2*L1*L2 +2*(L1+L2)+1 - max(L1,L2)   ) ! same as min(M1,M2)*(2*Lmin+1)*(Lmax+1)

write(*,'(A,I6)')               'Total number of molecular terms =   ', nTotal
write(*,*)

write(*,fmt='(A)',advance='no') 'Possible multiplicities are       = '
do i= 0, nint(2._dp*Smin)
 write(unit=*,fmt='(i3)',advance='no') nint( 2._dp*(Smax-Smin) ) + 2*i + 1   !nint( 2._dp*((Smax-Smin)+ real(i, dp) )+1._dp )
enddo
write(*,*)
write(*,*)

write(*,'(2A10,A12,2A4)') 'Lambda', 'n states', 'of which:', '+', '-'
do lambda=0, L1+L2
 nLambda = 0
  do l= abs(L1-L2), L1+L2
  do m=0, l
      if( m == lambda) nLambda = nLambda + 1
  enddo
  enddo
 write(*,'(2I10)',advance='no') lambda, nLambda
 if( lambda ==0) then
    nLambdaPlus  = int(nLambda/2)
    nLambdaMinus = int(nLambda/2)
    kk = p1*p2* (-1)**(L1+L2)
    if( nLambdaPlus+nLambdaMinus /= nLambda .and. kk==+1) nLambdaPlus  = nLambdaPlus+1
    if( nLambdaPlus+nLambdaMinus /= nLambda .and. kk==-1) nLambdaMinus = nLambdaMinus+1
    write(*,'(12x,2I4)',advance='no') nLambdaPlus, nLambdaMinus
 endif
 write(*,*)
enddo


end program WigmerWitmer

