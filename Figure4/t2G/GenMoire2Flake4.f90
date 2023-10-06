program GenMoire

   implicit none

   !double precision :: phideg = 21.051724d0 ! Rotation angle in degrees
   !double precision :: phideg = 6.01d0 ! Rotation angle in degrees
   !double precision :: phideg  = 0.5605482106766739d0, phideg2 = 1.12d0 ! Rotation angle in degrees
   !double precision :: phideg  =  -0.59d0, phideg2 = 1.13d0 ! Rotation angle in degrees
   double precision, parameter :: degtol12  =0.03d0, lattol12  = 0.002d0 ! Tolerance
   double precision, parameter :: degtol32 =0.03d0, lattol32 = 0.002d0 ! Tolerance
   double precision, parameter :: h = 3.35d0, dz = 35.0d0
  
   !!t2GBN 
   !double precision :: phideg  =  0.56d0, phideg2 = 1.13d0 ! Rotation angle in degrees
   !double precision, parameter :: aG = 2.46019d0, aBN0  = 2.505759d0 ! 2.46019d0 ! Lattice parameters
   !double precision, parameter :: aG2 = aG,       aBN02 = 2.46019d0 ! extep
   
!   !t3G
!   double precision :: phideg  =  1.12d0, phideg2 = 0.5605482106766739d0 ! Rotation angle in degrees
!   !double precision :: phideg  =  1.538d0, phideg2 = 1.1d0 ! Rotation angle in degrees
!   !double precision :: phideg  =  1.1d0, phideg2 = 1.1d0 ! Rotation angle in degrees
!   double precision, parameter :: aG = 2.46019d0, aBN0  = 2.46019d0 ! ! 2.46019d0 ! Lattice parameters
!   double precision, parameter :: aG2 = aG,       aBN02 = 2.46019d0 ! extep
   
   !t2BN- Circular Flake
   double precision :: phi12deg  =  2.0d0 !1.5d0 ! Substrate : Rotation angle in degrees
   double precision :: phi32deg  =  XXXd0 !1.5d0 ! Flake      : Rotation angle in degrees
   double precision :: Radius    =  20d0 !20d0 ! Flake Radius in aG
!!! G/G
!!   double precision, parameter :: aG = 2.46019d0, aBN0  = 2.46019d0 ! 2.46019d0 ! ! 2.46019d0 ! Lattice parameters
!!   character(*), parameter :: Sp_L1(2) = (/'C','C'/), Sp_L2(2) =(/'C','C'/)
!!   integer, parameter :: N_L1(2) = (/6,6/), N_L2(2) = (/6,6/)
!!! BN/BN
!!   double precision, parameter :: aG = 2.505759d0, aBN0  = 2.505759d0 ! 2.46019d0 ! ! 2.46019d0 ! Lattice parameters
!!   character(*), parameter :: Sp_L1(2) = (/'B','N'/), Sp_L2(2) =(/'B','N'/) 
!!   integer, parameter :: N_L1(2) = (/5,7/), N_L2(2) = (/5,7/)
!! G/BN
!   double precision, parameter :: aG = 2.505759d0 ! Lattice parameters :Reference layer
!   double precision, parameter :: aBN0  = 2.460190 ! Lattice parameters :Rotated layer
!   character(*), parameter :: Sp_L1(2) = (/'C','C'/), Sp_L2(2) =(/'B','N'/) ! L2 is Reference, Bottom Layer
!   integer, parameter :: N_L1(2) = (/6,6/), N_L2(2) = (/5,7/) ! L2 is Reference,Bottom Layer
! BN/BN
   !double precision, parameter :: aG = 2.505759d0 ! Lattice parameters :Reference Middle layer
   !double precision, parameter :: aBN3  = 2.505759d0 ! Lattice parameters :Rotated Top layer
   !double precision, parameter :: aBN1  = 2.505759d0 ! 2.460190 ! Lattice parameters :Rotated Bottom layer
   double precision, parameter :: aG = 2.46019d0 ! Lattice parameters :Reference Middle layer
   double precision, parameter :: aBN3  = 2.46019d0 ! Lattice parameters :Rotated Top layer
   double precision, parameter :: aBN1  = 2.46019d0 ! 2.460190 ! Lattice parameters :Rotated Bottom layer
   !character(*), parameter :: Sp_L1(2) = (/'B','N'/), Sp_L2(2) =(/'B','N'/), Sp_L3(2) =(/'C','C'/) ! L2 is Reference, Middle Layer
   !integer, parameter :: N_L1(2) = (/5,7/), N_L2(2) = (/5,7/), N_L3(2) = (/6,6/) ! L2 is Reference,Bottom Layer
   character(*), parameter :: Sp_L1(2) = (/'C','C'/), Sp_L2(2) =(/'C','C'/), Sp_L3(2) =(/'C','C'/) ! L2 is Reference, Middle Layer
   integer, parameter :: N_L1(2) = (/5,7/), N_L2(2) = (/5,7/), N_L3(2) = (/5,7/) ! L2 is Reference,Bottom Layer
   
   integer, parameter :: supkFinal = 2 !300 ! Even number 


   character(*), parameter :: select = 'delta'
   logical, parameter :: fdfbuild = .true. !.false. !.true.


   double precision, parameter :: basisA(2,2) = reshape((/ 0.0d0/3.0d0, 0.0d0/3.0d0, 1.0d0/3.0d0,1.0d0/3.0d0/),(/2,2/)) ! A position 
   double precision, parameter :: basisB(2,2) = reshape((/ 1.0d0/3.0d0, 1.0d0/3.0d0, 2.0d0/3.0d0,2.0d0/3.0d0/),(/2,2/)) ! B position
   double precision, parameter :: basisC(2,2) = reshape((/-1.0d0/3.0d0,-1.0d0/3.0d0, 0.0d0/3.0d0,0.0d0/3.0d0/),(/2,2/)) ! C position

   ! Define each layer basis using basisA, basisB, basisC
   double precision :: basis1(2,2) = basisA ! Rotated Bottom layer   !reshape((/0.0d0/3.0d0,0.0d0/3.0d0, 1.0d0/3.0d0,1.0d0/3.0d0/),(/2,2/)) ! Rotated Bottom layer
   double precision :: basis2(2,2) = basisA ! Reference layer        !reshape((/0.0d0/3.0d0,0.0d0/3.0d0, 1.0d0/3.0d0,1.0d0/3.0d0/),(/2,2/)) 
   double precision :: basis3(2,2) = basisA ! Rotated Flake layer    !reshape((/0.0d0/3.0d0,0.0d0/3.0d0, 1.0d0/3.0d0,1.0d0/3.0d0/),(/2,2/)) 

   double precision :: Rot(3,3) !, ucell0(2,2), ucellCirc(2,2), ucellF(2,2), basisCirc(2) ,bssCirc(2)

   double precision, parameter :: pi = acos(-1.0d0)

   double precision :: delta12, delta32, phi32, phi12, g12, g32, deltol12, deltol32, lambda1,lambda2,lambda3 !, delta, phi, g 
!   double precision :: deltabis, deltabis2, deltolbis, phibis, gbis, phibis2, gbis2
   integer :: n1, n2, nMoire1, nMoire2, nn12(4), n(6), ncell(3),m2(2,2),Rcell(2,2) !, nn0(4), nbis(4), n(4), ncell(2), nbis2(4)
   double precision, pointer :: X_L1(:,:), X_L2(:,:), X_L3(:,:)
   !double precision, pointer :: Xg1(:,:), Xbn1(:,:)
   integer :: i, supk !ncell,ncellbis, i !, a_rot1, b_rot1, a_ref, b_ref, a_rot2, b_rot2 

   double precision :: xShift, yShift, xx1, yy1, xx2, yy2, xx3, yy3
   double precision, parameter :: shift_to_frac(2,2)  = reshape((/1.0d0,0.0d0, -1.0d0,2.0d0/),(/2,2/))
   double precision :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, a1(3), a1Rot(3), a2(3), a2Rot(3), a3(3), a3Rot(3)

   integer :: firstN 

   n1 = 0
   n2 = 150
   nMoire1  = 0
   nMoire2  = 1
   delta12  = aBN1/aG - 1.0d0
   deltol12 = (aBN1 + lattol12)/aG - 1.0d0 - delta12
   call MoireFind(n1,n2,nMoire1,nMoire2,phi12deg,delta12,degtol12,deltol12,select,nn12)
!  32 10 32 11 32 10
!!22 21 21 22 22 21
!   nn0(1) = 22
!   nn0(2) = 21
!   nn0(3) = 21
!   nn0(4) = 22
!!17 16 16 17
!   nn0(1) = 17
!   nn0(2) = 16
!   nn0(3) = 16
!   nn0(4) = 17
!28 26 26 28
!30 29 31 28 33 27
!  nn0 = (/26,29,28,27/) !-2.4
!   nn0 = (/25, 29, 26, 28/)  ! -1.2
!   g = nn0(1)**2 + nn0(2)**2 + nn0(1)*nn0(2)
!   delta = sqrt(real(nn0(3)**2 + nn0(4)**2 + nn0(3)*nn0(4))/g)
!   phi = acos((2.0d0*nn0(1)*nn0(3)+2.0d0*nn0(2)*nn0(4) + nn0(1)*nn0(4) + nn0(2)*nn0(3))/(2.0d0*delta*g))
!
!   xx2    = 1.0d0
!   yy2    = 0.0d0
!   xx1    = (real(nn0(1)) * real(nn0(3)) + real(nn0(2))*real(nn0(4)) + real(nn0(2))*real(nn0(3)))/g
!   xx1    = xx1 + (real(nn0(1)) * real(nn0(4)) - real(nn0(2))*real(nn0(3)) )/g/2.d0
!   yy1    = (real(nn0(1)) * real(nn0(4)) - real(nn0(2))*real(nn0(3)) )/g*(sqrt(3.0d0)/2.0d0)
!
!   if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
!        phi = -phi
!   end if
!
!   write(*,'(a,f13.6)') '         Ang   :', phi*180.0d0/pi
!   write(*,'(a,f13.6)') '         Latt  :', aG*delta
!
!   write(*,*)
!   write(*,'(a,4i6)') "For 1st run, (a_rot1, b_rot1, a_ref, b_ref) is :", nn0  !nn0(1), nn0(2), nn0(3), nn0(4)
!   write(*,*)
!
!   deltabis  = aBN02/aG2 - 1.0d0
!   deltolbis = (aBN02 + lattol2)/aG2 - 1.0d0 - deltabis
!   supk = 50
!   call MoireFind2(n1,n2,nn0(3),nn0(4),supk,nMoire1,nMoire2,phideg2,deltabis,degtol2,deltolbis,select,nbis, supkFinal)
!
!   write(*,*)
!   write(*,'(a,4i6)') "For 2nd run, (a_rot2, b_rot2, a_ref, b_ref) is : ", nbis
!   write(*,'(a, i6)') "             supkFinal is                      : ", supkFinal
!   write(*,*)
!! 29    29    29    30    30    29   
!   n(1) = nn0(1)*supkFinal
!   n(2) = nn0(2)*supkFinal
!   n(3) = nn0(3)*supkFinal
!   n(4) = nn0(4)*supkFinal
!   n(5) = nbis(1)
!   n(6) = nbis(2)


   n(1:4)=nn12
   n(5:6)=n(3:4)
   n(1) = n(1)*supkFinal
   n(2) = n(2)*supkFinal
   n(3) = n(3)*supkFinal
   n(4) = n(4)*supkFinal
   n(5) = n(5)*supkFinal
   n(6) = n(6)*supkFinal

!basisCirc = basisCirc*n(3)*supkFinal

!basis2(1,1) = basis2(1,1) + real(supkFinal)/2.0d0
!basis2(2,1) = basis2(2,1) + real(supkFinal)/2.0d0
!basis2(1,2) = basis2(1,2) + real(supkFinal)/2.0d0
!basis2(2,2) = basis2(2,2) + real(supkFinal)/2.0d0

!basisCirc   = (/real(supkFinal)*0.5d0, real(supkFinal)*0.5d0/)
   write(*,*)
   write(*,'(a,2f17.9)') "basis1 ", basis1(1,1), basis1(2,1) !n(1), n(2), n(3), n(4), n(5), n(6)
   write(*,'(a,2f17.9)') "basis2 ", basis2(1,1), basis2(2,1) !n(1), n(2), n(3), n(4), n(5), n(6)
   write(*,'(a,2f17.9)') "basis3 ", basis3(1,1), basis3(2,1) !n(1), n(2), n(3), n(4), n(5), n(6)

!basisCirc(:) = basis2(:,1) 
!   write(*,'(a,2f17.9)') "basisF ", basisCirc(1), basisCirc(2) !n(1), n(2), n(3), n(4), n(5), n(6)

   write(*,*)
   write(*,'(a,6i6)') "Final : (a_rot1, b_rot1, a_ref, b_ref, a_rot3, b_rot3) is ", n !n(1), n(2), n(3), n(4), n(5), n(6)
   write(*,*)

   
      g12      = n(1)**2 + n(2)**2 + n(1)*n(2)
      delta12  = aBN1/aG ! sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
      phi12    = phi12deg*pi/180.0d0 ! acos((2.0d0*n(3)*n(1)+2.0d0*n(4)*n(2) + n(3)*n(2) + n(4)*n(1))/(2.0d0*delta*g))

      g32      = n(5)**2 + n(6)**2 + n(5)*n(6)
      delta32  = aBN3/aG ! sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
      phi32    = phi32deg*pi/180.0d0 ! acos((2.0d0*n(3)*n(1)+2.0d0*n(4)*n(2) + n(3)*n(2) + n(4)*n(1))/(2.0d0*delta*g))
!      gbis     = n(5)**2 + n(6)**2 + n(5)*n(6)
!      deltabis = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/gbis)
!      phibis   = acos((2.0d0*n(3)*n(5)+2.0d0*n(4)*n(6) + n(3)*n(6) + n(4)*n(5))/(2.0d0*deltabis*gbis))

      ncell(1) = g12
      ncell(2) = n(3)**2 + n(4)**2 + n(3)*n(4)
      ncell(3) = g32

!      xx2    = 1.0d0
!      yy2    = 0.0d0
!      xx1    = (real(n(1)) * real(n(3)) + real(n(2))*real(n(4)) + real(n(2))*real(n(3)))/g
!      xx1    = xx1 + (real(n(1)) * real(n(4)) - real(n(2))*real(n(3)) )/g/2.d0
!      yy1    = (real(n(1)) * real(n(4)) - real(n(2))*real(n(3)) )/g*(sqrt(3.0d0)/2.0d0)
!!      xx3    = (real(n(5)) * real(n(3)) + real(n(6))*real(n(4)) + real(n(2))*real(n(3)))/g
!!      xx3    = xx3 + (real(n(5)) * real(n(4)) - real(n(6))*real(n(3)) )/g/2.d0
!!      yy3    = (real(n(5)) * real(n(4)) - real(n(6))*real(n(3)) )/g*(sqrt(3.0d0)/2.0d0)
!
!      if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
!           phi = -phi
!      end if
!!      if (atan2(yy3,xx3) < atan2(yy2,xx2)) then
!!           phibis = -phibis
!!      end if
      write(*,*)
      write(*,'(a)')       'Final angles   :'
      write(*,'(a,f13.6)') '     L1        :', phi12*180.0d0/pi
      write(*,'(a,f13.6)') '     L2(ref)   :', 0.0d0*180.0d0/pi
      write(*,'(a,f13.6)') '     L3(Flake) :', phi32*180.0d0/pi
      write(*,*)
      write(*,'(a)')       'Final lattices :'
      write(*,'(a,2f13.6)') '    L1        :', aG  *delta12        
      write(*,'(a,2f13.6)') '    L2(ref)   :', aG          
      write(*,'(a,2f13.6)') '    L3(Flake) :', aG  *delta32

!      if (Sp_L1(1) .eq. 'B' .or. SP_L1(1) .eq. 'N') then
!          write(*,'(a,2f13.6)') '          L1   :', aG  *delta , 9.11d0*(aG*delta-2.505d0)*(aG*delta-2.505d0)*1.d3
!      else
!          write(*,'(a,2f13.6)') '          L1   :', aG  *delta , 11.44d0*(aG*delta-2.4602d0)*(aG*delta-2.4602d0)*1.d3
!      end if
!      if (Sp_L2(1) .eq. 'B' .or. SP_L2(1) .eq. 'N') then
!          write(*,'(a,2f13.6)') '          L2   :', aG   ,  9.11d0*(aG-2.505d0)*(aG-2.505d0)*1.d3
!      else
!          write(*,'(a,2f13.6)') '          L2   :', aG   , 11.44d0*(aG-2.4602d0)*(aG-2.4602d0)*1.d3
!      end if
!      if (Sp_L3(1) .eq. 'B' .or. SP_L3(1) .eq. 'N') then
!          write(*,'(a,2f13.6)') '          L3   :', aG2 *deltabis  ,  9.11d0*(aG2 *deltabis-2.505d0)*(aG2 *deltabis-2.505d0)*1.d3
!      else
!          write(*,'(a,2f13.6)') '          L3   :', aG2 *deltabis  , 11.44d0*(aG2 *deltabis-2.4602d0)*(aG2*deltabis-2.4602d0)*1.d3
!      end if
!      write(*,*)
!      write(*,'(a)')       'Final N cells  :'
!      write(*,'(a,i8)')    '     L1        :', ncell(1)
!      write(*,'(a,i8)')    '     L2(ref)   :', ncell(2)
!      write(*,'(a,i8)')    '     L3(Flake) :', ncell(3)
!      write(*,*)


   if (fdfbuild) then
      
      allocate(X_L1(2,ncell(1)*2))
      allocate(X_L2(2,ncell(2)*2))
      allocate(X_L3(2,ncell(3)*2))
      
!      ucellCirc(:,1) = (/cos(phi32),sin(phi32)/)
!      ucellCirc(:,2) = (/cos(phi32+pi/3.0d0),sin(phi32+pi/3)/)
!      ucellCirc = ucellCirc*delta32*aG
!
!      ucell0(:,1)    = (/1.0d0, 0.0d0/)
!      ucell0(:,2)    = (/cos(pi/3.0d0),sin(pi/3)/)
!      basisCirc = 0.5d0* ( (n(3)-n(4))*ucell0(:,1) + (n(3)+2*n(4))*ucell0(:,2) )
!      write(*,'(a,2f17.9)') "basisF ", basisCirc(1), basisCirc(2) !n(1), n(2), n(3), n(4), n(5), n(6)
!
!   m2(:,1) = (/n(3),n(4)/)
!   m2(:,2) = (/-n(4),n(3)+n(4)/)
!   Rcell = 0.0d0
!   Rcell(1:2,1:2) = matmul(ucell0,m2)
!
!
!   write(*,*)
!   write(*,'(a,2f16.9)')   'Center in Cart.  : ', aG*basisCirc(1), aG*basisCirc(2)
!      bssCirc(2) = Rcell(1,1)*(basisCirc(2)-basisCirc(1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
!      bssCirc(1) = (basisCirc(1)-bssCirc(2)*Rcell(1,2))/Rcell(1,1)
!      bssCirc = mod(bssCirc,1.0d0)
!      if (bssCirc(1)<0.0d0) bssCirc(1) = bssCirc(1) + 1.0d0
!      if (bssCirc(2)<0.0d0) bssCirc(2) = bssCirc(2) + 1.0d0
!   write(*,'(a,2f16.9)')   'Center in Frac.  : ', bssCirc(1), bssCirc(2)

!      ucellF(2,1) = Rcell(1,1)*(ucellCirc(2,1)-ucellCirc(1,1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
!      ucellF(1,1) = (ucellCirc(1,1)-ucellF(2,1)*Rcell(1,2))/Rcell(1,1)
!      ucellF(2,2) = Rcell(1,1)*(ucellCirc(2,2)-ucellCirc(1,2)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
!      ucellF(1,2) = (ucellCirc(1,2)-ucellF(2,2)*Rcell(1,2))/Rcell(1,1)

!      write(*,*)
!      write(*,'(a)')       'Final N atoms  :'
!      write(*,'(a,i8)')    '        Subs   :', ncell(2)*2
!      write(*,'(a,i8)')    '        Flake  :', ncell(1)*2
!      write(*,*)
   
 !     write(*,'(a,i8)')    '   supkFinal   :', supkFinal   
      write(*,*)
      write(*,'(a)')       'Constructing the layers'
      write(*,*)

      call MoireConstruct(     X_L1,n(1:2),basis1    , phi12 *180.0d0/pi ,delta12) ! Rotated Bottom layer 
      call MoireConstruct(     X_L2,n(3:4),basis2    , 0.0d0 *180.0d0/pi ,1.0d0)   ! Reference Middle layer
     
      ! Flake layer with no shift, no rotation, no delta correction 
      call MoireConstruct(     X_L3,n(5:6),basisA    , 0.0d0 *180.0d0/pi ,1.0d0)   ! Flake Top layer ; no shift, no rotation, no delta correction
    
      write(*,*)


      write(*,'(a,f13.6)')    '   Radius (Angstrom)   :', Radius*aG
      !call MoireWrite3L(Xg,Xbn,n,aG,-h,dz,'moire.fdf',grSp,bnSp,grN,bnN)
      !call MoireWrite4L(Xg,Xg1,Xbn,nbis,aG,h,dz,'moire.fdf',grSp,bnSp,grN,bnN)
      !call MoireWrite3Leach(X_L1, X_L2, X_L3, n, aG, h, dz, 'moire.fdf', Sp_L1, SP_L2, SP_L3, N_L1, N_L2, N_L3)
      !call MoireWrite2Leach(X_L1, X_L2, n, aG, h, dz, 'moire.fdf', Sp_L1, SP_L2, N_L1, N_L2)
      
      call MoireWrite3L_Hexa(X_L1, X_L2, X_L3, n, aG, h, dz, 'moire.fdf', Sp_L1, SP_L2, SP_L3, N_L1, N_L2, N_L3, &
                             & basis3-basisA, phi32, delta32, Radius) ! Additional Info for the flake 
                            ! (   sliding , twist angle, lattice mismatch , flake size )
      
      lambda1  = (delta12)*aG/sqrt(2.0d0*(delta12)*(1.0d0-cos(phi12))+(delta12-1.0d0)**2)
      !lambda2  = (delta32)*aG/sqrt(2.0d0*(delta32)*(1.0d0-cos(phi32))+(delta32-1.0d0)**2)
      write(*,'(a,3f13.6)') 'lambda from phi, delta :', lambda1 !, lambda2
      lambda1  = aG*(delta12*sqrt(real(ncell(1))))
      lambda2  = aG*sqrt(real(ncell(2)))
      !lambda3  = aG*delta32*sqrt(real(ncell(3)))
      write(*,'(a,3f13.6)') 'lambda from ncell      :', lambda1, lambda2 !, lambda3
      write(*,*)
      !write(*,'(i9)') ncell(1)*2+ncell(2)*2 
      !write(*,*)
      !write(*,*)
      write(*,'(a,3i9)') 'FixA: ', 1+ncell(1)*2+ncell(2)*2 !, 1+ncell(1)*2+ncell(2)*2
      write(*,'(a,i9,a,i9,a,i9,a,i9)') 'FixB: ', 1,',', ncell(1)*2,',', ncell(1)*2+1,',',ncell(1)*2 +ncell(2)*2 !, 1+ncell(1)*2+ncell(2)*2

!      write(*,'(a,3i9)') 'Fix: ', 90301, 1+ncell(1)*2 !, 1+ncell(1)*2+ncell(2)*2
   end if

contains

subroutine MoireFind(n1,n2,nM1,nM2,phiin,delta,degtol,lattol,select,mout)

   implicit none

   double precision :: pi = acos(-1.0d0)

   integer, intent(in) :: n1, n2, nM1, nM2
   double precision, intent(in) :: phiin, delta, degtol, lattol
   character(*), intent(in) :: select
   integer, intent(out) :: mout(4)

   integer :: a, b, ap, bp, m(2,2), mp(2,2), mm(2,2)
   integer :: nn, N, nn2
   double precision :: p, alpha, g, angletol, pcomp, d, phi
   double precision :: ep, ea,e, ep0, ea0, e0
   double precision :: xx1,yy1,xx2,yy2
   double precision :: m0p(4), m0a(4), m0(4)

   angletol = degtol*pi/180.0d0
   phi = phiin*pi/180.0d0
   ep0 = huge(1.0d0)
   ea0 = huge(1.0d0)
   e0 = huge(1.0d0)
   write(*,'(a,f16.9)') 'Angle:', phiin
   write(*,'(a,f16.9)') 'Delta:', delta
   !do ap=n1,n2; do bp=n1,ap; do a=n1,ap; do b=n1,ap
   do a=n1+1,n2; do b=n1,a; do ap=n1,n2; do bp=n1,n2
      m(:,1) = (/a,b/)
      m(:,2) = (/-b,a+b/)
      mp(:,1) = (/ap,bp/)
      mp(:,2) = (/-bp,ap+bp/)
      mm = m-mp
      !N = a**2 + b**2 + a*b      
      !if (mod(N,3)==0) write(*,*) N, a, b, ap, bp
      N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
      if (N>=nM1 .and. N<=nM2) then
         g = a**2 + b**2 + a*b
         p = sqrt(real(ap**2+bp**2+ap*bp)/g)
         nn = 2*a*ap + 2*b*bp + b*ap + a*bp
         alpha = acos(real(nn,8)/(2.0d0*p*g))

         xx2    = 1.0d0
         yy2    = 0.0d0
         xx1    = (real(a) * real(ap) + real(b)*real(bp) + real(b)*real(ap))/g
         xx1    = xx1 + (real(a) * real(bp) - real(b)*real(ap) )/g/2.d0
         yy1    = (real(a) * real(bp) - real(b)*real(ap) )/g*(sqrt(3.0d0)/2.0d0)

         if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
              alpha = -alpha
         end if

         !if (p < 1.0d0) then
         !   pcomp = 1.0d0/p
         !else
            pcomp = p
         !end if
         d = pcomp - 1.0d0
         ep = abs(d-delta)
         ea = abs(alpha-phi)
         if ((ep < lattol) .and. (ea < angletol)) then
         !if (abs(alpha-phi) < angletol) then
            if (ea < ea0) then
               ea0 = ea
               m0a = (/a,b,ap,bp/)
            end if
            if (ep < ep0) then
               ep0 = ep
               m0p = (/a,b,ap,bp/)
            end if
            !e = sqrt(ep**2/(delta**2+epsilon(1.0d0)) + ea**2/(phi**2+epsilon(1.0d0)))
            e = sqrt(ep**2/delta**2 + ea**2/phi**2)
            if (e < e0) then
               e0 = e
               m0 = (/a,b,ap,bp/)
            end if
            !write(*,'(4i6,2f16.9,i6)') a,b,ap,bp,pcomp, alpha*180.0d0/pi, N
         end if
      end if
   end do; end do; end do; end do
   write(*,*)
   write(*,'(a)') "                a     b     a'    b'       Delta            Angle       N   #1   #2"
   a = m0a(1)
   b = m0a(2)
   ap = m0a(3)
   bp = m0a(4)
   m(:,1) = (/a,b/)
   m(:,2) = (/-b,a+b/)
   mp(:,1) = (/ap,bp/)
   mp(:,2) = (/-bp,ap+bp/)
   mm = m-mp
   N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
   g = a**2 + b**2 + a*b
   p = sqrt(real(ap**2+bp**2+ap*bp)/g)
   nn = 2*a*ap + 2*b*bp + b*ap + a*bp
   alpha = acos(real(nn,8)/(2.0d0*p*g))
   xx2    = 1.0d0
   yy2    = 0.0d0
   xx1    = (real(a) * real(ap) + real(b)*real(bp) + real(b)*real(ap))/g
   xx1    = xx1 + (real(a) * real(bp) - real(b)*real(ap) )/g/2.d0
   yy1    = (real(a) * real(bp) - real(b)*real(ap) )/g*(sqrt(3.0d0)/2.0d0)
   if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
        alpha = -alpha
   end if
   !if (p < 1.0d0) then
   !   p = 1.0d0/p
   !   a = m0a(3)
   !   b = m0a(4)
   !   ap = m0a(1)
   !   bp = m0a(2)
   !end if
   d = p - 1.0d0
   nn = a**2 + a*b + b**2
   nn2 = ap**2 + ap*bp + bp**2
   write(*,'(a,4i6,2f16.9,i5,2i6)') 'For angle:  ',a,b,ap,bp,d,alpha*180.0d0/pi,N,nn, nn2
   a = m0p(1)
   b = m0p(2)
   ap = m0p(3)
   bp = m0p(4)
   m(:,1) = (/a,b/)
   m(:,2) = (/-b,a+b/)
   mp(:,1) = (/ap,bp/)
   mp(:,2) = (/-bp,ap+bp/)
   mm = m-mp
   N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
   g = a**2 + b**2 + a*b
   p = sqrt(real(ap**2+bp**2+ap*bp)/g)
   nn = 2*a*ap + 2*b*bp + b*ap + a*bp
   alpha = acos(real(nn,8)/(2.0d0*p*g))
   xx2    = 1.0d0
   yy2    = 0.0d0
   xx1    = (real(a) * real(ap) + real(b)*real(bp) + real(b)*real(ap))/g
   xx1    = xx1 + (real(a) * real(bp) - real(b)*real(ap) )/g/2.d0
   yy1    = (real(a) * real(bp) - real(b)*real(ap) )/g*(sqrt(3.0d0)/2.0d0)
   if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
        alpha = -alpha
   end if
   !if (p < 1.0d0) then
   !   p = 1.0d0/p
   !   a = m0p(3)
   !   b = m0p(4)
   !   ap = m0p(1)
   !   bp = m0p(2)
   !end if
   d = p - 1.0d0
   nn = a**2 + a*b + b**2
   nn2 = ap**2 + ap*bp + bp**2
   write(*,'(a,4i6,2f16.9,i5,2i6)') 'For delta:  ',a,b,ap,bp,d,alpha*180.0d0/pi,N, nn, nn2
   a = m0(1)
   b = m0(2)
   ap = m0(3)
   bp = m0(4)
   m(:,1) = (/a,b/)
   m(:,2) = (/-b,a+b/)
   mp(:,1) = (/ap,bp/)
   mp(:,2) = (/-bp,ap+bp/)
   mm = m-mp
   N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
   g = a**2 + b**2 + a*b
   p = sqrt(real(ap**2+bp**2+ap*bp)/g)
   nn = 2*a*ap + 2*b*bp + b*ap + a*bp
   alpha = acos(real(nn,8)/(2.0d0*p*g))
   xx2    = 1.0d0
   yy2    = 0.0d0
   xx1    = (real(a) * real(ap) + real(b)*real(bp) + real(b)*real(ap))/g
   xx1    = xx1 + (real(a) * real(bp) - real(b)*real(ap) )/g/2.d0
   yy1    = (real(a) * real(bp) - real(b)*real(ap) )/g*(sqrt(3.0d0)/2.0d0)
   if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
        alpha = -alpha
   end if
   !if (p < 1.0d0) then
   !   p = 1.0d0/p
   !   a = m0(3)
   !   b = m0(4)
   !   ap = m0(1)
   !   bp = m0(2)
   !end if
   d = p - 1.0d0
   nn = a**2 + a*b + b**2
   nn2 = ap**2 + ap*bp + bp**2
   write(*,'(a,4i6,2f16.9,i5,2i6)') 'For both:   ',a,b,ap,bp,d,alpha*180.0d0/pi,N,nn,nn2
   if (select=='angle' .or. select=='ANGLE') then
      mout = m0a
   else if (select=='delta' .or. select=='DELTA') then
      mout = m0p
   else if (select=='both' .or. select=='BOTH') then
      mout = m0
   end if
   if (all(mout==0)) stop 'All integers are 0'

end subroutine MoireFind

subroutine MoireFind2(n1,n2,n01,n02,supk,nM1,nM2,phiin,delta,degtol,lattol,select,mout, supkFinal)

   implicit none

   double precision :: pi = acos(-1.0d0)

   integer, intent(in) :: n1, n2, nM1, nM2, n01, n02, supk
   double precision, intent(in) :: phiin, delta, degtol, lattol
   character(*), intent(in) :: select
   integer, intent(out) :: mout(4), supkFinal
   integer :: supki,flag


   integer :: a, b, ap, bp, m(2,2), mp(2,2), mm(2,2)
   integer :: nn, N, nn2
   double precision :: p, alpha, g, angletol, pcomp, d, phi
   double precision :: ep, ea,e, ep0, ea0, e0
   double precision :: xx1,yy1,xx2,yy2 
   double precision :: m0p(4), m0a(4), m0(4)

   
   angletol = degtol*pi/180.0d0
   phi  = phiin*pi/180.0d0
   ep0  = huge(1.0d0)
   ea0  = huge(1.0d0)
   e0   = huge(1.0d0)
   write(*,'(a,f16.9)') 'Angle:', phiin
   write(*,'(a,f16.9)') 'Delta:', delta
   do supki=1,supk
      flag = 0
      ap = supki*n01
      bp = supki*n02
      print*, "supki=", supki
      !do a=n1,n2; do b=n1,a; do ap=n1,a; do bp=n1,a
      !do a=n1,ap-5; do b=n1,ap-5
      do a=n1,n2; do b=n1,n2
         !m(:,1) = (/a,b/)
         !m(:,2) = (/-b,a+b/)
         !mp(:,1) = (/ap,bp/)
         !mp(:,2) = (/-bp,ap+bp/)
         mp(:,1) = (/ap,bp/)
         mp(:,2) = (/-bp,ap+bp/)
         m(:,1) = (/a,b/)
         m(:,2) = (/-b,a+b/)
         mm = m-mp
         !N = a**2 + b**2 + a*b      
         !if (mod(N,3)==0) write(*,*) N, a, b, ap, bp
         N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
         if (N>=nM1 .and. N<=nM2) then
            g = a**2 + b**2 + a*b
            p = sqrt(real(ap**2+bp**2+ap*bp)/g)
            nn = 2*a*ap + 2*b*bp + b*ap + a*bp
            alpha = acos(real(nn,8)/(2.0d0*p*g))
   xx2    = 1.0d0
   yy2    = 0.0d0
   xx1    = (real(a) * real(ap) + real(b)*real(bp) + real(b)*real(ap))/g
   xx1    = xx1 + (real(a) * real(bp) - real(b)*real(ap) )/g/2.d0
   yy1    = (real(a) * real(bp) - real(b)*real(ap) )/g*(sqrt(3.0d0)/2.0d0)
   if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
        alpha = -alpha
   end if
            !if (p < 1.0d0) then
            !   pcomp = 1.0d0/p
            !else
               pcomp = p
            !end if
            d = pcomp - 1.0d0
            ep = abs(d-delta)
            ea = abs(alpha-phi)
            if ((ep < lattol) .and. (ea < angletol)) then
               flag = 1
            !if (abs(alpha-phi) < angletol) then
               if (ea < ea0) then
                  ea0 = ea
                  m0a = (/a,b,ap,bp/)
               end if
               if (ep < ep0) then
                  ep0 = ep
                  m0p = (/a,b,ap,bp/)
               end if
               !e = sqrt(ep**2/(delta**2+epsilon(1.0d0)) + ea**2/(phi**2+epsilon(1.0d0)))
               e = sqrt(ep**2/delta**2 + ea**2/phi**2)
               if (e < e0) then
                  e0 = e
                  m0 = (/a,b,ap,bp/)
               end if
               !write(*,'(4i6,2f16.9,i6)') a,b,ap,bp,pcomp, alpha*180.0d0/pi, N
            end if
         end if
      end do; end do!; end do; end do
      if (flag==1) then
          supkFinal = supki
          exit
      end if
   end do
   write(*,*)
   write(*,'(a)') "                a     b     a'    b'       Delta            Angle       N   #1   #2"
   a = m0a(1)
   b = m0a(2)
   ap = m0a(3)
   bp = m0a(4)
   m(:,1) = (/a,b/)
   m(:,2) = (/-b,a+b/)
   mp(:,1) = (/ap,bp/)
   mp(:,2) = (/-bp,ap+bp/)
   mm = m-mp
   N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
   g = a**2 + b**2 + a*b
   p = sqrt(real(ap**2+bp**2+ap*bp)/g)
   nn = 2*a*ap + 2*b*bp + b*ap + a*bp
   alpha = acos(real(nn,8)/(2.0d0*p*g))
   xx2    = 1.0d0
   yy2    = 0.0d0
   xx1    = (real(a) * real(ap) + real(b)*real(bp) + real(b)*real(ap))/g
   xx1    = xx1 + (real(a) * real(bp) - real(b)*real(ap) )/g/2.d0
   yy1    = (real(a) * real(bp) - real(b)*real(ap) )/g*(sqrt(3.0d0)/2.0d0)
   if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
        alpha = -alpha
   end if
   !if (p < 1.0d0) then
   !   p = 1.0d0/p
   !   a = m0a(3)
   !   b = m0a(4)
   !   ap = m0a(1)
   !   bp = m0a(2)
   !end if
   d = p - 1.0d0
   nn = a**2 + a*b + b**2
   nn2 = ap**2 + ap*bp + bp**2
   write(*,'(a,4i6,2f16.9,i5,2i6)') 'For angle:  ',a,b,ap,bp,d,alpha*180.0d0/pi,N,nn, nn2
   a = m0p(1)
   b = m0p(2)
   ap = m0p(3)
   bp = m0p(4)
   m(:,1) = (/a,b/)
   m(:,2) = (/-b,a+b/)
   mp(:,1) = (/ap,bp/)
   mp(:,2) = (/-bp,ap+bp/)
   mm = m-mp
   N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
   g = a**2 + b**2 + a*b
   p = sqrt(real(ap**2+bp**2+ap*bp)/g)
   nn = 2*a*ap + 2*b*bp + b*ap + a*bp
   alpha = acos(real(nn,8)/(2.0d0*p*g))
   xx2    = 1.0d0
   yy2    = 0.0d0
   xx1    = (real(a) * real(ap) + real(b)*real(bp) + real(b)*real(ap))/g
   xx1    = xx1 + (real(a) * real(bp) - real(b)*real(ap) )/g/2.d0
   yy1    = (real(a) * real(bp) - real(b)*real(ap) )/g*(sqrt(3.0d0)/2.0d0)
   if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
        alpha = -alpha
   end if
   !if (p < 1.0d0) then
   !   p = 1.0d0/p
   !   a = m0p(3)
   !   b = m0p(4)
   !   ap = m0p(1)
   !   bp = m0p(2)
   !end if
   d = p - 1.0d0
   nn = a**2 + a*b + b**2
   nn2 = ap**2 + ap*bp + bp**2
   write(*,'(a,4i6,2f16.9,i5,2i6)') 'For delta:  ',a,b,ap,bp,d,alpha*180.0d0/pi,N, nn, nn2
   a = m0(1)
   b = m0(2)
   ap = m0(3)
   bp = m0(4)
   m(:,1) = (/a,b/)
   m(:,2) = (/-b,a+b/)
   mp(:,1) = (/ap,bp/)
   mp(:,2) = (/-bp,ap+bp/)
   mm = m-mp
   N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
   g = a**2 + b**2 + a*b
   p = sqrt(real(ap**2+bp**2+ap*bp)/g)
   nn = 2*a*ap + 2*b*bp + b*ap + a*bp
   alpha = acos(real(nn,8)/(2.0d0*p*g))
   xx2    = 1.0d0
   yy2    = 0.0d0
   xx1    = (real(a) * real(ap) + real(b)*real(bp) + real(b)*real(ap))/g
   xx1    = xx1 + (real(a) * real(bp) - real(b)*real(ap) )/g/2.d0
   yy1    = (real(a) * real(bp) - real(b)*real(ap) )/g*(sqrt(3.0d0)/2.0d0)
   if (atan2(yy1,xx1) < atan2(yy2,xx2)) then
        alpha = -alpha
   end if
   !if (p < 1.0d0) then
   !   p = 1.0d0/p
   !   a = m0(3)
   !   b = m0(4)
   !   ap = m0(1)
   !   bp = m0(2)
   !end if
   d = p - 1.0d0
   nn = a**2 + a*b + b**2
   nn2 = ap**2 + ap*bp + bp**2
   write(*,'(a,4i6,2f16.9,i5,2i6)') 'For both:   ',a,b,ap,bp,d,alpha*180.0d0/pi,N,nn,nn2
   if (select=='angle' .or. select=='ANGLE') then
      mout = m0a
   else if (select=='delta' .or. select=='DELTA') then
      mout = m0p
   else if (select=='both' .or. select=='BOTH') then
      mout = m0
   end if
   if (all(mout==0)) stop 'All integers are 0'

end subroutine MoireFind2

subroutine MoireConstruct(X,mm,bss,angle,p)

   implicit none

   double precision, parameter :: pi = acos(-1.0d0)

   double precision, intent(out) :: X(:,:)
   integer, intent(in) :: mm(2)
   double precision, intent(in) :: bss(:,:), angle, p

   integer :: ncell, nbasis, i1, i2, j, k, n1, n2, mint(2,2)
   double precision :: ucell(2,2), Rcell(2,2), a, m(2,2)
   double precision :: f(2), sq
   double precision, pointer :: basis(:,:)
   integer, save :: nlayer=0

   print*, "constructing layer with angle = ", angle

   ncell = mm(1)**2 + mm(1)*mm(2) + mm(2)**2
   nbasis = size(bss,2)
   allocate(basis(2,nbasis))
   if (size(X) /= 2*ncell*nbasis) then
      stop 'Wrong size of X array'
   end if
   a = angle*pi/180.0d0
   ucell(:,1) = (/cos(a),sin(a)/)
   ucell(:,2) = (/cos(a+pi/3.0d0),sin(a+pi/3)/)
   ucell = ucell*p
   m(:,1) = (/mm(1),mm(2)/)
   m(:,2) = (/-mm(2),mm(1)+mm(2)/)
   Rcell = matmul(ucell,m)
   basis = matmul(ucell,bss)
   sq = sqrt(real(ncell))
   mint = m
   n1 = gcd(mint(1,1),mint(2,1))
   n2 = ncell/n1
   nlayer = nlayer  + 1
   write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'Size of layer ',nlayer,':   ',n1,' x ',n2, '   (',mm(1),' x ',mm(2),')'
   !write(*,*) ucell
   k = 0
   !write(*,*) sqrt(Rcell(1,1)**2+Rcell(2,1)**2)/sqrt(ucell(1,1)**2+ucell(2,1)**2)
   do i1 =1,n1; do i2=1,n2
      do j=1,nbasis
         k = k + 1
         X(:,k) = (i1-1)*ucell(:,1) + (i2-1)*ucell(:,2) + basis(:,j)
      end do
   end do; end do
   do i1=1,ncell*nbasis
      !print*, Rcell(1,1), X(2,i1), X(1,i1)
      f(2) = Rcell(1,1)*(X(2,i1)-X(1,i1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      f(1) = (X(1,i1)-f(2)*Rcell(1,2))/Rcell(1,1)
      X(:,i1) = mod(f,1.0d0)
      if (X(1,i1)<0.0d0) X(1,i1) = X(1,i1) + 1.0d0
      if (X(2,i1)<0.0d0) X(2,i1) = X(2,i1) + 1.0d0
   end do
   deallocate(basis)

end subroutine MoireConstruct

subroutine MoireConstructFlake(X,mm,bss,angle,p)

   implicit none

   double precision, parameter :: pi = acos(-1.0d0)

   double precision, intent(out) :: X(:,:)
   integer, intent(in) :: mm(2) 
   double precision, intent(in) :: bss(:,:), angle, p !, bssC(:,:) 

   integer :: ncell, nbasis, i1, i2, j, k, n1, n2, mint(2,2)
   double precision :: ucell(2,2), ucell0(2,2),a, m(2,2), Rcell0(3,3) , Rcell(2,2) !, a, m(2,2) 
   double precision :: f(2), sq, dist
   double precision, pointer :: basis(:,:)
   integer, save :: nlayer=0

   print*, "constructing layer with angle = ", angle

   ncell =  ( mm(1)**2 + mm(1)*mm(2) + mm(2)**2 )
   nbasis = size(bss,2)
   allocate(basis(2,nbasis))
   if (size(X) /= 2*ncell*nbasis) then
      stop 'Wrong size of X array'
   end if
   a = angle*pi/180.0d0
   ucell(:,1) = (/cos(a),sin(a)/)
   ucell(:,2) = (/cos(a+pi/3.0d0),sin(a+pi/3)/)
   ucell = ucell*p
   basis = matmul(ucell,bss)

   ucell0(:,1) = (/1.0d0,0.0d0/)
   ucell0(:,2) = (/cos(pi/3.0d0),sin(pi/3)/)
   m(:,1) = (/mm(1),mm(2)/)
   m(:,2) = (/-mm(2),mm(1)+mm(2)/)
   Rcell = matmul(ucell0,m)

!   a = angle*pi/180.0d0
!   ucell(:,1) = (/cos(a),sin(a)/)
!   ucell(:,2) = (/cos(a+pi/3.0d0),sin(a+pi/3)/)
!   !ucell = ucell*p
!   basis = matmul(ucell,bss)
!
!   Rcell = matmul(ucell,m)

   !write(*,*) bssC(1), bssC(2)
   sq = sqrt(real(ncell))
   mint = m
   n1 = gcd(mint(1,1),mint(2,1))
   n2 = ncell/n1
   nlayer = nlayer  + 1
   write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'Size of layer ',nlayer,':   ',n1,' x ',n2, '   (',mm(1),' x ',mm(2),')'
   !write(*,*) ucell
   k = 0
   !write(*,*) sqrt(Rcell(1,1)**2+Rcell(2,1)**2)/sqrt(ucell(1,1)**2+ucell(2,1)**2)
   do i1 =1,n1; do i2=1,n2
      do j=1,nbasis
         k = k + 1
         X(:,k) = (i1-1)*ucell(:,1) + (i2-1)*ucell(:,2)  + basis(:,j) !&
         f  = (i1-1)*ucell(:,1) + (i2-1)*ucell(:,2)  + basis(:,j)  & !Rcell(1,1)*(Xs3(1,i))+Rcell(1,2)*(Xs3(2,i)) -basisCirc(1)
              & - 0.5d0*(mm(1)-mm(2))*ucell(:,1) - 0.5d0*(mm(1)+2.d0*mm(2))*ucell(:,2) 
         dist  = sqrt( f(1)*f(1) + f(2)*f(2) )
         !if (dist <= RRadius) then
         !k = k + 1
         !X(:,k) = (i1-1)*ucell(:,1) + (i2-1)*ucell(:,2)  + basis(:,j) !&
         
         !end if
!                  & - 0.5d0*(mm(1)-mm(2))*ucell(:,1) - 0.5d0*(mm(1)+2.d0*mm(2))*ucell(:,2)  !(0.5d0)*Rcell(:,1) + (0.5d0)*Rcell(:,2)
         !XX      = (i1-1)*ucell(1,1) + (i2-1)*ucell(1,2)  + basis(1,j) &
         !         & - 0.5d0*(mm(1)-mm(2))*ucell(:,1) - 0.5d0*(mm(1)+2.d0*mm(2))*ucell(:,2)  !(0.5
      end do
   end do; end do
!   a1 = Rcell(:,1) ! (/Rcell(1,1), Rcell(2,1)/)
!   a2 = Rcell(:,2) !(/Rcell(1,2), Rcell(2,2),0.0d0/) !Rcell(:,2)
!   !a3 = (/Rcell(1,1), Rcell(2,1),/) !z
!   !Rot = reshape((/cos(-a),sin(-a),0.0d0,-sin(-a),cos(-a),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
!   Rot = reshape((/cos(-a),sin(-a),-sin(-a),cos(-a)/),(/2,2/))
!   a1Rot = matmul(Rot,a1)
!   a2Rot = matmul(Rot,a2)
 

   do i1=1,ncell*nbasis
      !print*, Rcell(1,1), X(2,i1), X(1,i1)
      f(2) = Rcell(1,1)*(X(2,i1)-X(1,i1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      f(1) = (X(1,i1)-f(2)*Rcell(1,2))/Rcell(1,1)
      !f(2) = a1Rot(1)*(X(2,i1)-X(1,i1)*a1Rot(2)/a1Rot(1))/(a2Rot(2)*a1Rot(1)-a2Rot(1)*a1Rot(2))
      !f(1) = (X(1,i1)-f(2)*a2Rot(1))/a1Rot(1)
      !X(:,i1) = f !mod(f,1.0d0)
      X(:,i1) = mod(f,1.0d0)
      if (X(1,i1)<0.0d0) X(1,i1) = X(1,i1) + 1.0d0
      if (X(2,i1)<0.0d0) X(2,i1) = X(2,i1) + 1.0d0
      X(1,i1) = X(1,i1)*p
      X(2,i1) = X(2,i1)*p
   end do
   deallocate(basis)

end subroutine MoireConstructFlake

subroutine MoireConstructShift(X,mm,bss,angle,p,yShift)

   implicit none

   double precision, parameter :: pi = acos(-1.0d0)

   double precision, intent(out) :: X(:,:)
   integer, intent(in) :: mm(2)
   double precision, intent(in) :: bss(:,:), angle, p

   integer :: ncell, nbasis, i1, i2, j, k, n1, n2, mint(2,2)
   double precision :: ucell(2,2), Rcell(2,2), a, m(2,2)
   double precision :: f(2), sq
   double precision, pointer :: basis(:,:)
   integer, save :: nlayer=0

   double precision :: yShift

   ncell = mm(1)**2 + mm(1)*mm(2) + mm(2)**2
   nbasis = size(bss,2)
   allocate(basis(2,nbasis))
   if (size(X) /= 2*ncell*nbasis) then
      stop 'Wrong size of X array'
   end if
   a = angle*pi/180.0d0
   ucell(:,1) = (/cos(a),sin(a)/)
   ucell(:,2) = (/cos(a+pi/3.0d0),sin(a+pi/3)/)
   ucell = ucell*p
   m(:,1) = (/mm(1),mm(2)/)
   m(:,2) = (/-mm(2),mm(1)+mm(2)/)
   Rcell = matmul(ucell,m)
   basis = matmul(ucell,bss)
   sq = sqrt(real(ncell))
   mint = m
   n1 = gcd(mint(1,1),mint(2,1))
   n2 = ncell/n1
   nlayer = nlayer  + 1
   write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'Size of layer ',nlayer,':   ',n1,' x ',n2, '   (',mm(1),' x ',mm(2),')'
   !write(*,*) ucell
   k = 0
   !write(*,*) sqrt(Rcell(1,1)**2+Rcell(2,1)**2)/sqrt(ucell(1,1)**2+ucell(2,1)**2)
   do i1 =1,n1; do i2=1,n2
      do j=1,nbasis
         k = k + 1
         X(:,k) = (i1-1)*ucell(:,1) + (i2-1)*ucell(:,2) + basis(:,j)
      end do
   end do; end do
   print*, yShift
   do i1=1,ncell*nbasis
      !print*, Rcell(1,1), X(2,i1), X(1,i1)
      f(2) = Rcell(1,1)*(X(2,i1)-X(1,i1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      f(1) = (X(1,i1)-f(2)*Rcell(1,2))/Rcell(1,1)
      !X(:,i1) = mod(f,1.0d0)
      !if (X(1,i1)<0.0d0) X(1,i1) = X(1,i1) + 1.0d0
      !if (X(2,i1)<0.0d0) X(2,i1) = X(2,i1) + 1.0d0

      X(:,i1) = [mod(f(1),1.0), mod(f(2),1.0)+yShift]!, 0.5_dp+d]
      if (X(1,i1)<0.0) X(1,i1) = mod(X(1,i1) + 1.0,1.0)
      if (X(2,i1)<0.0) X(2,i1) = mod(X(2,i1) + 1.0,1.0)
      if (X(1,i1)>1.0) X(1,i1) = mod(X(1,i1) - 1.0,1.0)
      if (X(2,i1)>1.0) X(2,i1) = mod(X(2,i1) - 1.0,1.0)

   end do
   deallocate(basis)

end subroutine MoireConstructShift

subroutine MoireWrite(Xs,Xo,mm,a,h,z,name,sp1,sp2,N1,N2)

   implicit none

   integer, parameter :: u=10, u2=11, u3=12, u4=13, u5=14, u6=15

   double precision, intent(in) :: Xs(:,:), Xo(:,:), a, h, z
   integer, intent(in) :: mm(4), N1(:), N2(:)
   character(*), intent(in) :: name, sp1(:), sp2(:)

   double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2)
   double precision :: rz
   integer :: sz1, sz2, i, j, k, nb1, nb2
   character(60), pointer :: Sp(:)
   integer, pointer :: Nel(:), num1(:), num2(:)
   double precision :: phirad, e1(3)

   open(u,FILE=name,STATUS='replace')
   open(u2,FILE="BLBL.frac",STATUS='replace')
   open(u3,FILE="BLBL.basis1",STATUS='replace')
   open(u4,FILE="BLBL.basis2",STATUS='replace')
   open(u5,FILE="BLBL.mol", STATUS='replace')
   open(u6,FILE="sublattices.dat", STATUS='replace')
   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)
   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
   Rcell = 0.0d0
   Rcell(1:2,1:2) = matmul(ucell,m1)
   Rcell(3,3) = z
   sz1 = size(Xs,2)
   sz2 = size(Xo,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
   allocate(Sp(nb1+nb2))
   allocate(Nel(nb1+nb2))
   allocate(num1(nb1))
   allocate(num2(nb2))
   Sp(1) = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = 1
out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
   end do out2
   write(u,'(a,i0)') 'NumberOfSpecies', k
   write(u,'(a)') '%block ChemicalSpeciesLabel'
   do i=1,k
      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
   end do
   !a1 = (/Rcell(1,1), Rcell(1,2), Rcell(1,3)/)
   !a2 = (/Rcell(2,1), Rcell(2,2), Rcell(2,3)/)
   !a3 = (/Rcell(3,1), Rcell(3,2), Rcell(3,3)/)
   a1 = Rcell(:,1)
   a2 = Rcell(:,2)
   a3 = Rcell(:,3)
   e1 = (/1,0,0/)
   phirad = -acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
   Rot = reshape((/cos(phirad),sin(phirad),0.0d0,-sin(phirad),cos(phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
   a1Rot = matmul(Rot,a1)
   a2Rot = matmul(Rot,a2)
   a3Rot = matmul(Rot,a3)
   !print*, a1, a2, a3
   !print*, a1Rot, a2Rot, a3Rot
   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
   write(u,'(a,i0)') 'NumberOfAtoms  ', sz1+sz2
   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz2)
   write(u,'(a)') 'LatticeConstant  1.00 Ang'
   write(u,'(a)') '%block LatticeVectors'
   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
   write(u,'(a)') '%endblock LatticeVectors'
   write(u,*)
   write(u2,*)
   xlo = 0.0d0
   xhi = a1Rot(1)
   ylo = 0.0d0
   yhi = a2Rot(2)
   zlo = 0.0d0
   zhi = 35.0d0
   xy = a2Rot(1) - 0.000000001d0
   xz = 0.0d0
   yz = 0.0d0
   write(u2,'(a,i0)') 'Lattice vectors:'
   write(u2,'(a,3f16.9)') 'a:', xhi, 0.0d0, 0.0d0
   write(u2,'(a,3f16.9)') 'b:', xy, yhi, 0.0d0
   write(u2,'(a,3f16.9)') 'c:', 0.0d0, 0.0d0, 35.0d0
   !write(u2,'(a,i0)') 'Lattice vectors:'
   !do i=1,3
   !   write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   !end do
   write(u2,*)
   write(u2,'(a)') 'Fractional coordinates:'
   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   write(u5,*)
   write(u5,'(i8,a)') sz1+sz2, ' atoms'
   write(u5,'(i1,a)') 2, ' atom types'
   write(u5,*)
   write(u5,'(2f16.9,a)') xlo, xhi, ' xlo xhi'
   write(u5,'(2f16.9,a)') ylo, yhi, ' ylo yhi'
   write(u5,'(2f16.9,a)') zlo, zhi, ' zlo zhi'
   write(u5,'(3f16.9,a)') xy, xz, yz, ' xy xz yz'
   write(u5,*)
   write(u5,'(a)') ' Masses'
   write(u5,*)
   write(u5,'(i1, 3f16.9)') 1, 12.0107
   write(u5,'(i1, 3f16.9)') 2, 12.0107
   write(u5,*)
   write(u5,'(a)') ' Atoms'
   write(u5,*)
   rz = (z/2.0d0 - h/2.0d0)/z
   j = 0
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xs(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', 2, ' 0.0 ', Xs(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   rz = (z/2.0d0 + h/2.0d0)/z
   do i=1,sz2
      j = j+1
      write(u,'(3f16.9,i8)') Xo(:,i), rz, num2(mod(i+1,nb2)+1)
      write(u2,'(A8,3f16.9)') 'C', Xo(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xo(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', 1, ' 0.0 ', Xo(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   close(u)

end subroutine MoireWrite

subroutine MoireWrite3L(Xs,Xo,mm,a,h,z,name,sp1,sp2,N1,N2)

   implicit none

   integer, parameter :: u=10, u2=11, u3=12, u4=13, u5=14, u6=15

   double precision, intent(in) :: Xs(:,:), Xo(:,:), a, h, z
   integer, intent(in) :: mm(4), N1(:), N2(:)
   character(*), intent(in) :: name, sp1(:), sp2(:)

   double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2)
   double precision :: rz
   integer :: sz1, sz2, i, j, k, nb1, nb2
   character(60), pointer :: Sp(:)
   integer, pointer :: Nel(:), num1(:), num2(:)
   double precision :: phirad, e1(3)

   open(u,FILE=name,STATUS='replace')
   open(u2,FILE="BLBL.frac",STATUS='replace')
   open(u3,FILE="BLBL.basis1",STATUS='replace')
   open(u4,FILE="BLBL.basis2",STATUS='replace')
   open(u5,FILE="BLBL.mol", STATUS='replace')
   open(u6,FILE="sublattices.dat", STATUS='replace')
   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)
   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
   Rcell = 0.0d0
   Rcell(1:2,1:2) = matmul(ucell,m1)
   Rcell(3,3) = z
   sz1 = size(Xs,2)
   sz2 = size(Xo,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
   allocate(Sp(nb1+nb2))
   allocate(Nel(nb1+nb2))
   allocate(num1(nb1))
   allocate(num2(nb2))
   Sp(1) = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = 1
out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
   end do out2
   write(u,'(a,i0)') 'NumberOfSpecies', k
   write(u,'(a)') '%block ChemicalSpeciesLabel'
   do i=1,k
      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
   end do
   !a1 = (/Rcell(1,1), Rcell(1,2), Rcell(1,3)/)
   !a2 = (/Rcell(2,1), Rcell(2,2), Rcell(2,3)/)
   !a3 = (/Rcell(3,1), Rcell(3,2), Rcell(3,3)/)
   a1 = Rcell(:,1)
   a2 = Rcell(:,2)
   a3 = Rcell(:,3)
   e1 = (/1,0,0/)
   phirad = -acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
   Rot = reshape((/cos(phirad),sin(phirad),0.0d0,-sin(phirad),cos(phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
   a1Rot = matmul(Rot,a1)
   a2Rot = matmul(Rot,a2)
   a3Rot = matmul(Rot,a3)
   !print*, a1, a2, a3
   !print*, a1Rot, a2Rot, a3Rot
   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
   write(u,'(a,i0)') 'NumberOfAtoms  ', sz1+sz2+sz1
   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz2+sz1)
   write(u,'(a)') 'LatticeConstant  1.00 Ang'
   write(u,'(a)') '%block LatticeVectors'
   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
   write(u,'(a)') '%endblock LatticeVectors'
   write(u,*)
   write(u2,*)
   xlo = 0.0d0
   xhi = a1Rot(1)
   ylo = 0.0d0
   yhi = a2Rot(2)
   zlo = 0.0d0
   zhi = 35.0d0
   xy = a2Rot(1)
   xz = 0.0d0
   yz = 0.0d0
   write(u2,'(a,i0)') 'Lattice vectors:'
   write(u2,'(a,3f16.9)') 'a:', xhi, 0.0d0, 0.0d0
   write(u2,'(a,3f16.9)') 'b:', xy, yhi, 0.0d0
   write(u2,'(a,3f16.9)') 'c:', 0.0d0, 0.0d0, 35.0d0
   !write(u2,'(a,i0)') 'Lattice vectors:'
   !do i=1,3
   !   write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   !end do
   write(u2,*)
   write(u2,'(a)') 'Fractional coordinates:'
   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   write(u5,*)
   write(u5,'(i8,a)') sz1+sz2+sz1, ' atoms'
   write(u5,'(i1,a)') 2, ' atom types'
   write(u5,*)
   write(u5,'(2f16.9,a)') xlo, xhi, ' xlo xhi'
   write(u5,'(2f16.9,a)') ylo, yhi, ' ylo yhi'
   write(u5,'(2f16.9,a)') zlo, zhi, ' zlo zhi'
   write(u5,'(3f16.9,a)') xy, xz, yz, ' xy xz yz'
   write(u5,*)
   write(u5,'(a)') ' Masses'
   write(u5,*)
   write(u5,'(i1, 3f16.9)') 1, 12.0107
   write(u5,'(i1, 3f16.9)') 2, 12.0107
   write(u5,*)
   write(u5,'(a)') ' Atoms'
   write(u5,*)
   rz = (z/2.0d0 - h/2.0d0)/z
   j = 0
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xs(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', 2, ' 0.0 ', Xs(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   rz = (z/2.0d0 + h/2.0d0)/z
   do i=1,sz2
      j = j+1
      write(u,'(3f16.9,i8)') Xo(:,i), rz, num2(mod(i+1,nb2)+1)
      write(u2,'(A8,3f16.9)') 'C', Xo(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xo(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', 1, ' 0.0 ', Xo(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   rz = (z/2.0d0 + 3.0*h/2.0d0)/z
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xs(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', 2, ' 0.0 ', Xs(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   close(u)

end subroutine MoireWrite3L

subroutine MoireWrite2Leach(Xs1,Xs2,mm,a,h,z,fdfname,sp1,sp2,N1,N2)

   implicit none

   integer, parameter :: u=10, u2=11, u3=12, u4=13, u5=14, u6=15

   double precision, intent(in) :: Xs1(:,:), Xs2(:,:),  a, h, z
   integer, intent(in) :: mm(4), N1(:), N2(:)
   character(*), intent(in) :: fdfname, sp1(:), sp2(:)

   double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2), m3(2,2)
   double precision :: rz
   integer :: sz1, sz2, i, j, k, nb1, nb2 !, atomSp(3)=(/5,7,6/) !corresponding (/'B','N','C'/)
   character(60), pointer :: Sp(:)
   double precision, pointer :: SpMass(:)
   character :: atomSp(3)=(/'B','N','C'/)
   double precision :: atomMass(3)=(/10.800000191, 14.000000000, 12.010000229/)
   integer, pointer :: Nel(:),  num1(:), num2(:) 
   double precision :: phirad,phirad2, e1(3)

   open(u,FILE=fdfname,STATUS='replace')
   open(u2,FILE="BLBL.frac",STATUS='replace')
   open(u3,FILE="BLBL.basis1",STATUS='replace')
   open(u4,FILE="BLBL.basis2",STATUS='replace')
   open(u5,FILE="BLBL.mol", STATUS='replace')
   open(u6,FILE="sublattices.dat", STATUS='replace')

   

   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)

   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
!   m3(:,1) = (/mm(5),mm(6)/)
!   m3(:,2) = (/-mm(6),mm(5)+mm(6)/)
   Rcell = 0.0d0
   Rcell(1:2,1:2) = matmul(ucell,m2)
   Rcell(3,3) = z
   sz1 = size(Xs1,2)
   sz2 = size(Xs2,2)
!   sz3 = size(Xs3,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
!   nb3 = size(sp3)

   allocate(Sp(nb1+nb2))
   allocate(Nel(nb1+nb2))
   allocate(SpMass(nb1+nb2))
   allocate(num1(nb1))
   allocate(num2(nb2))
!   allocate(num3(nb3))

   Sp(1)  = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = 1
   if (Sp(1)==atomSp(1)) SpMass(1)=atomMass(1)
   if (Sp(1)==atomSp(2)) SpMass(1)=atomMass(2)
   if (Sp(1)==atomSp(3)) SpMass(1)=atomMass(3)
   print*, SpMass(k)

out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      print*, SpMass(k)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
   end do out2
!out3:do i=1,nb3
!      do j=1,k
!         if (sp3(i)==Sp(j)) then
!            num3(i) = j
!            cycle out3
!         end if
!      end do
!      k = k+1
!      num3(i) = k
!      Sp(k) = sp3(i)
!      Nel(k) = N3(i)
!      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
!      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
!      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
!   end do out3

   write(u,'(a,i0)') 'NumberOfSpecies', k
   write(u,'(a)') '%block ChemicalSpeciesLabel'
   do i=1,k
      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
   end do
   !a1 = (/Rcell(1,1), Rcell(1,2), Rcell(1,3)/)
   !a2 = (/Rcell(2,1), Rcell(2,2), Rcell(2,3)/)
   !a3 = (/Rcell(3,1), Rcell(3,2), Rcell(3,3)/)
   a1 = Rcell(:,1)
   a2 = Rcell(:,2)
   a3 = Rcell(:,3)
   e1 = (/1,0,0/)
   phirad = acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
   write(*,'(f16.9, a)')  phirad*180.0d0/pi , ' : Angle of commensurate cell'
   Rot = reshape((/cos(-phirad),sin(-phirad),0.0d0,-sin(-phirad),cos(-phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
   a1Rot = matmul(Rot,a1)
   a2Rot = matmul(Rot,a2)
   a3Rot = matmul(Rot,a3)
   !print*, a1, a2, a3
   !print*, a1Rot, a2Rot, a3Rot
   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
   write(u,'(a,i0)') 'NumberOfAtoms  ', (sz1+sz2)
   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz2)
   write(u,'(a)') 'LatticeConstant  1.00 Ang'
   write(u,'(a)') '%block LatticeVectors'
   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
   write(u,'(a)') '%endblock LatticeVectors'
   write(u,*)
   write(u2,*)
   xlo = 0.0d0
   xhi = a1Rot(1)
   ylo = 0.0d0
   yhi = a2Rot(2)
   zlo = 0.0d0
   zhi = 35.0d0
   xy = a2Rot(1) - 0.000000001d0
   xz = 0.0d0
   yz = 0.0d0
   write(*,'(3f16.9)') (a1Rot(j),j=1,2)
   write(*,'(3f16.9)') (a2Rot(j),j=1,2)
   write(u2,'(a,i0)') 'Lattice vectors:'
   write(u2,'(a,3f16.9)') 'a:', xhi, 0.0d0, 0.0d0
   write(u2,'(a,3f16.9)') 'b:', xy, yhi, 0.0d0
   write(u2,'(a,3f16.9)') 'c:', 0.0d0, 0.0d0, 35.0d0
   !write(u2,'(a,i0)') 'Lattice vectors:'
   !do i=1,3
   !   write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   !end do
   write(u2,*)
   write(u2,'(a)') 'Fractional coordinates:'
   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   write(u5,*)
   write(u5,'(i8,a)') (sz1+sz2), ' atoms'
   write(u5,'(i1,a)') k, ' atom types'
   write(u5,*)
   write(u5,'(2f16.9,a)') xlo, xhi, ' xlo xhi'
   write(u5,'(2f16.9,a)') ylo, yhi, ' ylo yhi'
   write(u5,'(2f16.9,a)') zlo, zhi, ' zlo zhi'
   write(u5,'(3f16.9,a)') xy, xz, yz, ' xy xz yz'
   write(u5,*)
   write(u5,'(a)') ' Masses'
   write(u5,*)
   do j=1,k
   write(u5,'(i1, 3f16.9)') j, SpMass(j)
   write(*,'(i1, 3f16.9)') j, SpMass(j)
   end do   
   write(u5,*)
   write(u5,'(a)') ' Atoms'
   write(u5,*)
   rz = (z/2.0d0 - 1.0d0*h/2.0d0)/z
   j = 0
   write(*,'(f16.3)') rz*z
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs1(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') sp1(mod(i+1,nb1)+1), Xs1(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs1(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', num1(mod(i+1,nb1)+1), ' 0.0 ', Xs1(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8,i8)') j, 1, mod(i+1,nb1)+1
   end do
   rz = (z/2.0d0 + 1.0d0*h/2.0d0)/z
   write(*,'(f16.3)') rz*z
   do i=1,sz2
      j = j+1
      write(u,'(3f16.9,i8)') Xs2(:,i), rz, num2(mod(i+1,nb2)+1)
      write(u2,'(A8,3f16.9)') sp2(mod(i+1,nb2)+1), Xs2(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs2(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', num2(mod(i+1,nb2)+1), ' 0.0 ', Xs2(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8,i8)') j, 2, mod(i+1,nb2)+1+2
   end do
!   rz = (z/2.0d0 + 1.0d0*h/2.0d0)/z
!   write(*,'(f16.3)') rz*z
!   do i=1,sz3
!      j = j+1
!      write(u,'(3f16.9,i8)') Xs3(:,i), rz, num3(mod(i+1,nb3)+1)
!      write(u2,'(A8,3f16.9)') sp3(mod(i+1,nb3)+1), Xs3(:,i), rz
!      write(u3,'(a,3f16.9,a)') 'basis ', Xs3(:,i), rz, ' &'
!      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
!      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
!      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', num3(mod(i+1,nb3)+1), ' 0.0 ', Xs3(:,i), rz, 0, 0, 0 
!      write(u6,'(i5,i8,i8)') j, 3, mod(i+1,nb3)+1+4
!   end do
   write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   close(u)

end subroutine MoireWrite2Leach

subroutine MoireWrite3L_Circ(Xs1,Xs2,Xs3,mm,a,h,z,fdfname,sp1,sp2,sp3,N1,N2,N3, phi32, delta32, rradius )

   implicit none

   !integer, parameter :: u=10, u2=11, u3=12, u4=13, u5=14, u6=15, u7=16
   integer, parameter :: u5=14, u6=15, u7=16

   double precision, intent(in) :: Xs1(:,:), Xs2(:,:), Xs3(:,:),  a, h, z, rradius, phi32, delta32 !, basisCirc(:) !, ucellF(:,:)
   integer, intent(in) :: mm(6), N1(:), N2(:), N3(:)
   character(*), intent(in) :: fdfname, sp1(:), sp2(:), sp3(:)

   double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2) , m3(2,2), Rcell0(3,3)
   double precision :: a1(3),a2(3),a3(3),e1(3), a1Rot(3),a2Rot(3),a3Rot(3) 
   double precision :: rz,rz1,rz2, ang
   integer :: sz1, sz2, sz3, i, j, k, j2, nb1, nb2,nb3, nn1,nn2,nn !, atomSp(3)=(/5,7,6/) !corresponding (/'B','N','C'/)
   character(60), pointer :: Sp(:)
   double precision, pointer :: SpMass(:)
   character :: atomSp(4)=(/'H','B','N','C'/)
   double precision :: atomMass(4)=(/1.00784, 10.800000191, 14.000000000, 12.010000229/)
   integer, pointer :: Nel(:),  num1(:), num2(:) , num3(:)
   double precision :: phirad, dist, vecs(3), f(2), basisCirc(2) !phirad2, e1(3), dist, vecs(3), f(2) , basisCirc(2)

!   open(u,FILE=fdfname,STATUS='replace')
!   open(u2,FILE="BLBL.frac",STATUS='replace')
!   open(u3,FILE="BLBL.basis1",STATUS='replace')
!   open(u4,FILE="BLBL.basis2",STATUS='replace')
   open(u5,FILE="BLBL.mol", STATUS='replace')
   open(u6,FILE="sublattices.dat", STATUS='replace')
   open(u7,FILE="layerIndex.dat", STATUS='replace')

   

!   ucell(:,1) = (/1.0d0,0.0d0/)
!   ucell(:,2) = (/0.5d0,0.5d0*sqrt(3.0d0)/)

   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)

   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
   m3(:,1) = (/mm(5),mm(6)/)
   m3(:,2) = (/-mm(6),mm(5)+mm(6)/)
   Rcell0= 0.0d0
   Rcell = 0.0d0
   Rcell0(1:2,1:2) = matmul(ucell,m2)
   Rcell0(3,3) = z

   a1 = Rcell0(:,1)
   a2 = Rcell0(:,2)
   a3 = Rcell0(:,3)
   e1 = (/1,0,0/)
   phirad = acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
   write(*,'(f16.9, a)')  phirad*180.0d0/pi , ' : Angle of commensurate cell'
   Rot = reshape((/cos(-phirad),sin(-phirad),0.0d0,-sin(-phirad),cos(-phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
   Rcell(:,1) = matmul(Rot,a1)
   Rcell(:,2) = matmul(Rot,a2)
   Rcell(:,3) = matmul(Rot,a3)
   a1Rot = matmul(Rot,a1)
   a2Rot = matmul(Rot,a2)
   a3Rot = matmul(Rot,a3)

   write(*,*)
   basisCirc = 0.0d0 
   basisCirc(:) = 0.5d0*Rcell(1:2,1)+0.5d0*Rcell(1:2,2)
   write(*,'(a,2f16.9)')   'Center in Cart.  : ', basisCirc(1), basisCirc(2)
!      bssCirc(2) = Rcell(1,1)*(basisCirc(2)-basisCirc(1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
!      bssCirc(1) = (basisCirc(1)-bssCirc(2)*Rcell(1,2))/Rcell(1,1)
! !     bssCirc = mod(bssCirc,1.0d0)
! !     if (bssCirc(1)<0.0d0) bssCirc(1) = bssCirc(1) + 1.0d0
! !     if (bssCirc(2)<0.0d0) bssCirc(2) = bssCirc(2) + 1.0d0
!  write(*,'(a,2f16.9)')   'Center in Frac.  : ', bssCirc(1), bssCirc(2)



!   write(*,'(a,2f16.9)')   'Center in Frac.  : ', bssCirc(1), bssCirc(2)
!   basisCirc(1) = Rcell(1,1)*bssCirc(1)+Rcell(1,2)*bssCirc(2)
!   basisCirc(2) = Rcell(2,1)*bssCirc(1)+Rcell(2,2)*bssCirc(2)
!   write(*,'(a,2f16.9)')   'Center in Cart.  : ', basisCirc(1), basisCirc(2)


      !f(2) = Rcell(1,1)*(X(2,i1)-X(1,i1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      !f(1) = (X(1,i1)-f(2)*Rcell(1,2))/Rcell(1,1)
      !X(:,i1) = mod(f,1.0d0)
      !if (X(1,i1)<0.0d0) X(1,i1) = X(1,i1) + 1.0d0
      !if (X(2,i1)<0.0d0) X(2,i1) = X(2,i1) + 1.0d0

!      ucellF(2,1) = Rcell(1,1)*(ucellCirc(2,1)-ucellCirc(1,1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
!      ucellF(1,1) = (ucellCirc(1,1)-ucellF(2,1)*Rcell(1,2))/Rcell(1,1)
!      ucellF(:,1) = mod(ucellF(:,1),1.0d0)
!      if (ucellF(1,1)<0.0d0) ucellF(1,1) = ucellF(1,1) + 1.0d0
!      if (ucellF(2,1)<0.0d0) ucellF(2,1) = ucellF(2,1) + 1.0d0
 !     ucellF(2,2) = Rcell(1,1)*(ucellCirc(2,2)-ucellCirc(1,2)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
 !     ucellF(1,2) = (ucellCirc(1,2)-ucellF(2,2)*Rcell(1,2))/Rcell(1,1)
 !     ucellF(:,2) = mod(ucellF(:,2),1.0d0)
 !     if (ucellF(1,2)<0.0d0) ucellF(1,2) = ucellF(1,2) + 1.0d0
 !     if (ucellF(2,2)<0.0d0) ucellF(2,2) = ucellF(2,2) + 1.0d0

   sz1 = size(Xs1,2)
   sz2 = size(Xs2,2)
   sz3 = size(Xs3,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
   nb3 = size(sp3)

   allocate(Sp(nb1+nb2+nb3))
   allocate(Nel(nb1+nb2+nb3))
   allocate(SpMass(nb1+nb2+nb3))
   allocate(num1(nb1))
   allocate(num2(nb2))
   allocate(num3(nb3))

!Sp(1) = atomSp(1)
!Nel(1)= 1
!SpMass(1)=atomMass(1)

!   Sp(2)  = sp1(1)
!   Nel(2) = N1(1)
!   k = 2
   Sp(1)  = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = k
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      if (Sp(k)==atomSp(4)) SpMass(k)=atomMass(4)
   print*, SpMass(k)

out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      if (Sp(k)==atomSp(4)) SpMass(k)=atomMass(4)
      print*, SpMass(k)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      if (Sp(k)==atomSp(4)) SpMass(k)=atomMass(4)
   end do out2
out3:do i=1,nb3
      do j=1,k
         if (sp3(i)==Sp(j)) then
            num3(i) = j
            cycle out3
         end if
      end do
      k = k+1
      num3(i) = k
      Sp(k) = sp3(i)
      Nel(k) = N3(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      if (Sp(k)==atomSp(4)) SpMass(k)=atomMass(4)
   end do out3


!   write(u,'(a,i0)') 'NumberOfSpecies', k
!   write(u,'(a)') '%block ChemicalSpeciesLabel'
!   do i=1,k
!      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
!   end do
   !a1 = (/Rcell(1,1), Rcell(1,2), Rcell(1,3)/)
   !a2 = (/Rcell(2,1), Rcell(2,2), Rcell(2,3)/)
   !a3 = (/Rcell(3,1), Rcell(3,2), Rcell(3,3)/)
!   a1 = Rcell(:,1)
!   a2 = Rcell(:,2)
!   a3 = Rcell(:,3)
!   e1 = (/1,0,0/)
!   phirad = acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
!   write(*,'(f16.9, a)')  phirad*180.0d0/pi , ' : Angle of commensurate cell'
!   Rot = reshape((/cos(-phirad),sin(-phirad),0.0d0,-sin(-phirad),cos(-phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
!   a1Rot = matmul(Rot,a1)
!   a2Rot = matmul(Rot,a2)
!   a3Rot = matmul(Rot,a3)
   
!   write(*,'(a,2f16.9)')   'Center in Frac.  : ', bssCirc(1), bssCirc(2)
!   !basisCirc(1) = a1Rot(1)*bssCirc(1)+a2Rot(1)*bssCirc(2)
!   !basisCirc(2) = a1Rot(2)*bssCirc(1)+a2Rot(2)*bssCirc(2)
!   basisCirc(1) = Rcell(1,1)*bssCirc(1)+Rcell(1,2)*bssCirc(2)
!   basisCirc(2) = Rcell(2,1)*bssCirc(1)+Rcell(2,2)*bssCirc(2)
!   write(*,'(a,2f16.9)')   'Center in Cart.  : ', basisCirc(1), basisCirc(2)
   
!   ucellF(2,1) = a1Rot(1)*(ucellCirc(2,1)-ucellCirc(1,1)*a1Rot(2)/a1Rot(1))/(a2Rot(2)*a1Rot(1)-a2Rot(1)*a1Rot(2))
!   ucellF(1,1) = (ucellCirc(1,1)-ucellF(2,1)*a2Rot(1))/a1Rot(1)
!   ucellF(2,2) = a1Rot(1)*(ucellCirc(2,2)-ucellCirc(1,2)*a1Rot(2)/a1Rot(1))/(a2Rot(2)*a1Rot(1)-a2Rot(1)*a1Rot(2))
!   ucellF(1,2) = (ucellCirc(1,2)-ucellF(2,2)*a2Rot(1))/a1Rot(1)

   !print*, a1, a2, a3
   !print*, a1Rot, a2Rot, a3Rot
!   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
!   write(u,'(a,i0)') 'NumberOfAtoms  ', (sz1+sz2+sz3)
!   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz2+sz3)
!   write(u,'(a)') 'LatticeConstant  1.00 Ang'
!   write(u,'(a)') '%block LatticeVectors'
!   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
!   write(u,'(a)') '%endblock LatticeVectors'
!   write(u,*)
!   write(u2,*)
   xlo = 0.0d0
   xhi = a1Rot(1)
   ylo = 0.0d0
   yhi = a2Rot(2)
   zlo = 0.0d0
   zhi = 35.0d0
   xy = a2Rot(1) - 0.000000001d0
   xz = 0.0d0
   yz = 0.0d0
   write(*,'(3f16.9)') (a1Rot(j),j=1,2)
   write(*,'(3f16.9)') (a2Rot(j),j=1,2)
!   write(u2,'(a,i0)') 'Lattice vectors:'
!   write(u2,'(a,3f16.9)') 'a:', xhi, 0.0d0, 0.0d0
!   write(u2,'(a,3f16.9)') 'b:', xy, yhi, 0.0d0
!   write(u2,'(a,3f16.9)') 'c:', 0.0d0, 0.0d0, 35.0d0
   !write(u2,'(a,i0)') 'Lattice vectors:'
   !do i=1,3
   !   write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   !end do
!   write(u2,*)
!   write(u2,'(a)') 'Fractional coordinates:'
!   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
!   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   write(u5,*)
   write(u5,'(a)') 'NNN atoms'
   write(u5,'(i1,a)') k, ' atom types'
   write(u5,*)
   write(u5,'(2f16.9,a)') xlo, xhi, ' xlo xhi'
   write(u5,'(2f16.9,a)') ylo, yhi, ' ylo yhi'
   write(u5,'(2f16.9,a)') zlo, zhi, ' zlo zhi'
   write(u5,'(3f16.9,a)') xy, xz, yz, ' xy xz yz'
   write(u5,*)
   write(u5,'(a)') ' Masses'
   write(u5,*)
   do j=1,k
   write(u5,'(i1, 3f16.9)') j, SpMass(j)
   write(*,'(i1, 3f16.9)') j, SpMass(j)
   end do   
   write(u5,*)
   write(u5,'(a)') ' Atoms'
   write(u5,*)
   j = 0
   rz = (z/2.0d0 - 3.0d0*h/2.0d0)/z
   write(*,'(f16.3)') rz*z
   do i=1,sz1
      j = j+1
      !write(u,'(3f16.9,i8)') Xs1(:,i), rz, num1(mod(i+1,nb1)+1)
      !write(u2,'(A8,3f16.9)') sp1(mod(i+1,nb1)+1), Xs1(:,i), rz
      vecs(1)  = Rcell(1,1)*(Xs1(1,i))+Rcell(1,2)*(Xs1(2,i)) 
      vecs(2)  = Rcell(2,1)*(Xs1(1,i))+Rcell(2,2)*(Xs1(2,i)) 
      vecs(3)  = rz*z
      !write(u3,'(a,3f16.9,a)') 'basis ', Xs1(:,i), rz, ' &'
      !write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', num1(mod(i+1,nb1)+1), ' 0.0 ', vecs(:), 0,0,0 !Xs1(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j,  mod(i+1,nb1)+1
      write(u7,'(i5,i8)') j,  1 
   end do
   rz1 = rz
      write(*,*)
      write(*,'(a)')       'Final N atoms  :'
      write(*,'(a,i8)')    '     L1        :', j
   j2 = 0
   rz = (z/2.0d0 - 1.0d0*h/2.0d0)/z
!   write(*,'(f16.3)') rz*z
   do i=1,sz2
      j = j +1
      j2= j2+1
      !write(u,'(3f16.9,i8)') Xs2(:,i), rz, num2(mod(i+1,nb2)+1)
      !write(u2,'(A8,3f16.9)') sp2(mod(i+1,nb2)+1), Xs2(:,i), rz
      vecs(1)  = Rcell(1,1)*(Xs2(1,i))+Rcell(1,2)*(Xs2(2,i)) 
      vecs(2)  = Rcell(2,1)*(Xs2(1,i))+Rcell(2,2)*(Xs2(2,i)) 
      vecs(3)  = rz*z
      !write(u3,'(a,3f16.9,a)') 'basis ', Xs2(:,i), rz, ' &'
      !write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', num2(mod(i+1,nb2)+1), ' 0.0 ', vecs(:), 0,0,0 ! Xs2(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j,  mod(i+1,nb2)+1+2
      write(u7,'(i5,i8)') j,  2 
   end do
   rz2 = rz
      write(*,'(a,i8)')    '     L2(ref)   :', j2
   j2 = 0
   rz = (z/2.0d0 + 1.0d0*h/2.0d0)/z
 do i=1,sz3
      f(1)  = delta32 *( Rcell(1,1)*(Xs3(1,i))+Rcell(1,2)*(Xs3(2,i)) -basisCirc(1) )
      f(2)  = delta32 *( Rcell(2,1)*(Xs3(1,i))+Rcell(2,2)*(Xs3(2,i)) -basisCirc(2) )
      dist     = sqrt( f(1)*f(1) + f(2)*f(2) )
      if (dist <= rradius*a) then
         j = j+1
         j2 = j2+1
         
         vecs(1) = cos(phi32)*f(1) - sin(phi32)*f(2) + basisCirc(1)
         vecs(2) = sin(phi32)*f(1) + cos(phi32)*f(2) + basisCirc(2)
         vecs(3) = rz*z
         !write(u,'(3f16.9,i8)') Xs3(:,i)+bssCirc(:), rz, num3(mod(i+1,nb3)+1)
         !write(u,'(3f16.9,i8)') Xs3(:,i), rz, num3(mod(i+1,nb3)+1)
         !write(u2,'(A8,3f16.9)') sp3(mod(i+1,nb3)+1), Xs3(:,i), rz
         !write(u3,'(a,3f16.9,a)') 'basis ', Xs3(:,i), rz, ' &'
         !write(u4,'(a,i5,a,a)') 'basis ', j, ' 3', ' &'
         !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
         write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', num3(mod(i+1,nb3)+1), ' 0.0 ', vecs(:), 0,0,0 ! Xs3(:,i), rz, 0, 0, 0 
         !write(u6,'(i5,i8,i8)') j, 3, mod(i+1,nb3)+1
         write(u6,'(i5,i8)') j,  mod(i+1,nb3)+1+4
         write(u7,'(i5,i8)') j,  3 

!      else
!         if (mod(i+1,nb3)+1 .eq. 1) then ! A sublattice
!         do nn = 1,3 ! For the three nearest B sublattices
!            if (nn.eq.1) then 
!               nn1= 0.0d0; nn2= 0.0d0  
!            else if (nn.eq.2) then 
!               nn1=-1.0d0; nn2= 0.0d0 
!            else if (nn.eq.3) then 
!               nn1= 0.0d0; nn2=-1.0d0  
!            end if
!            vecs2(1)  = -basisCirc(1) + Rcell(1,1)*(Xs3(1,i+1))+Rcell(1,2)*(Xs3(2,i+1)) +nn1*ucellF(1,1)+nn2*ucellF(1,2) !-basisCirc(1)
!            vecs2(2)  = -basisCirc(2) + Rcell(2,1)*(Xs3(1,i+1))+Rcell(2,2)*(Xs3(2,i+1)) +nn1*ucellF(2,1)+nn2*ucellF(2,2) !-basisCirc(2)
!            dist2    = sqrt( vecs2(1)*vecs2(1) + vecs2(2)*vecs2(2) )
!            if (dist2 <= rradius*a) then
!               j = j+1
!               j2 = j2+1
!               write(u,'(3f16.9,i8)') Xs3(:,i), rz,1  !num1(mod(i+1,nb1)+1)
!               write(u2,'(A8,3f16.9)') '  H ' , Xs3(:,i), rz
!               write(u3,'(a,3f16.9,a)') 'basis ', Xs3(:,i), rz, ' &'
!               write(u4,'(a,i5,a,a)') 'basis ', j, ' 3', ' &'
!               !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
!               write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 3 ', 1, ' 0.0 ', Xs3(:,i), rz, 0, 0, 0 
!               write(u6,'(i5,i8,i8)') j, 3, mod(i+1,nb3)+1
!            end if
!         end do
!         else ! B sublattice
!         do nn = 1,3 ! For the three nearest A sublattices
!            if (nn.eq.1) then 
!               nn1= 0.0d0; nn2= 0.0d0  
!            else if (nn.eq.2) then 
!               nn1= 1.0d0; nn2= 0.0d0 
!            else if (nn.eq.3) then 
!               nn1= 0.0d0; nn2= 1.0d0  
!            end if
!            vecs2(1)  = -basisCirc(1) + Rcell(1,1)*(Xs3(1,i-1))+Rcell(1,2)*(Xs3(2,i-1)) +nn1*ucellF(1,1)+nn2*ucellF(1,2) !-basisCirc(1)
!            vecs2(2)  = -basisCirc(2) + Rcell(2,1)*(Xs3(1,i-1))+Rcell(2,2)*(Xs3(2,i-1)) +nn1*ucellF(2,1)+nn2*ucellF(2,2) !-basisCirc(2)
!            dist2    = sqrt( vecs2(1)*vecs2(1) + vecs2(2)*vecs2(2) )
!            if (dist2 <= rradius*a) then
!               j = j+1
!               j2 = j2+1
!               write(u,'(3f16.9,i8)') Xs3(:,i), rz, 1!num1(mod(i+1,nb1)+1)
!               write(u2,'(A8,3f16.9)') '  H ' , Xs3(:,i), rz
!               write(u3,'(a,3f16.9,a)') 'basis ', Xs3(:,i), rz, ' &'
!               write(u4,'(a,i5,a,a)') 'basis ', j, ' 3', ' &'
!               !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
!               write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 3 ', 1, ' 0.0 ',  Xs3(:,i), rz, 0, 0, 0 
!               write(u6,'(i5,i8,i8)') j, 3, mod(i+1,nb3)+1
!            end if
!         end do
!         end if
      end if
   end do
      write(*,'(a,i8)')    '     L3(Flake) :', j2
      write(*,'(a,i8)')    '        Total  :', j
      write(*,*)
   write(*,'(f16.3)') rz1*z
   write(*,'(f16.3)') rz2*z
   write(*,'(f16.3)') rz*z
      write(*,*)
   !write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   !close(u)
   !close(u2)
   close(u5)
   close(u6)
   close(u7)

end subroutine MoireWrite3L_Circ

subroutine MoireWrite3L_Hexa(Xs1,Xs2,Xs3,mm,a,h,z,fdfname,sp1,sp2,sp3,N1,N2,N3, basis32, phi32, delta32, rradius )

   implicit none

   !integer, parameter :: u=10, u2=11, u3=12, u4=13, u5=14, u6=15, u7=16
   integer, parameter :: u5=14, u6=15, u7=16

   double precision, intent(in) :: Xs1(:,:), Xs2(:,:), Xs3(:,:),  a, h, z, rradius, phi32, delta32, basis32(:,:) !, basisCirc(:) !, ucellF(:,:)
   integer, intent(in) :: mm(6), N1(:), N2(:), N3(:)
   character(*), intent(in) :: fdfname, sp1(:), sp2(:), sp3(:)

   !double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2) , m3(2,2), Rcell0(3,3), ucell0(2,2)
   double precision :: ucell(3,3), Rcell(3,3), m1(2,2), m2(2,2) , m3(2,2), Rcell0(3,3), ucell0(3,3)
   double precision :: a1(3),a2(3),a3(3),e1(3), a1Rot(3),a2Rot(3),a3Rot(3) 
   double precision :: rz,rz1,rz2, ang
   integer :: sz1, sz2, sz3, i, j, k, j2,j3, nb1, nb2,nb3, nn1,nn2,nn !, atomSp(3)=(/5,7,6/) !corresponding (/'B','N','C'/)
   character(60), pointer :: Sp(:)
   double precision, pointer :: SpMass(:)
   character :: atomSp(4)=(/'H','B','N','C'/)
   double precision :: atomMass(4)=(/1.00784, 10.800000191, 14.000000000, 12.010000229/)
   integer, pointer :: Nel(:),  num1(:), num2(:) , num3(:)
   double precision :: phirad, dist, vecs(3), f(2), basisCirc(2), basisCirc0(2) !,xx,yy  !phirad2, e1(3), dist, vecs(3), f(2) , basisCirc(2)

!   open(u,FILE=fdfname,STATUS='replace')
!   open(u2,FILE="BLBL.frac",STATUS='replace')
!   open(u3,FILE="BLBL.basis1",STATUS='replace')
!   open(u4,FILE="BLBL.basis2",STATUS='replace')
   open(u5,FILE="BLBL.mol", STATUS='replace')
   open(u6,FILE="sublattices.dat", STATUS='replace')
   open(u7,FILE="layerIndex.dat", STATUS='replace')

  
   ! ucell0 : reference layer lattice vectors 
   ucell0      = 0.0d0
   ucell0(1:2,1) = (/a,0.0d0/)
   ucell0(1:2,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)
   ucell0(3,3) = z ! Just for the dimension

   ! Rcell0 : Commensurate lattice vectors, which has the same orientation of ucell0
   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
   m3(:,1) = (/mm(5),mm(6)/)
   m3(:,2) = (/-mm(6),mm(5)+mm(6)/)
   Rcell0= 0.0d0
   Rcell0(1:2,1:2) = matmul(ucell0(1:2,1:2),m2)
   Rcell0(3,3) = z

   ! Rcell : Commensurate lattice vectors, which are aligned to x axis
   a1 = Rcell0(:,1)
   a2 = Rcell0(:,2)
   a3 = Rcell0(:,3)
   e1 = (/1,0,0/)
   phirad = acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
   write(*,'(f16.9, a)')  phirad*180.0d0/pi , ' : Angle of commensurate cell'
   Rot = reshape((/cos(-phirad),sin(-phirad),0.0d0,-sin(-phirad),cos(-phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
   Rcell(:,1) = matmul(Rot,a1)
   Rcell(:,2) = matmul(Rot,a2)
   Rcell(:,3) = matmul(Rot,a3)
   a1Rot = matmul(Rot,a1)
   a2Rot = matmul(Rot,a2)
   a3Rot = matmul(Rot,a3)

   ! ucell : reference layer lattice vectors, which has the same orientation of Rcell
   a1 = ucell0(:,1)
   a2 = ucell0(:,2)
   a3 = ucell0(:,3)
   ucell(:,1) = matmul(Rot,a1)
   ucell(:,2) = matmul(Rot,a2)
   ucell(:,3) = matmul(Rot,a3)

   ! basisCirc : rotation center in cartisian coordinate
   write(*,*)
   basisCirc = 0.0d0 
   basisCirc(:)  = 0.5d0*Rcell(1:2,1)+0.5d0*Rcell(1:2,2)
   basisCirc0(:) = 0.5d0*Rcell0(1:2,1)+0.5d0*Rcell0(1:2,2)
   write(*,'(a,2f16.9)')   'Center in Cart.  : ', basisCirc(1), basisCirc(2)

   sz1 = size(Xs1,2)
   sz2 = size(Xs2,2)
   sz3 = size(Xs3,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
   nb3 = size(sp3)

   allocate(Sp(nb1+nb2+nb3))
   allocate(Nel(nb1+nb2+nb3))
   allocate(SpMass(nb1+nb2+nb3))
   allocate(num1(nb1))
   allocate(num2(nb2))
   allocate(num3(nb3))

! To make it Hydrogen as the first atom species
Sp(1)    = atomSp(1)
Nel(1)   = 1
SpMass(1)= atomMass(1)

Sp(2)    = sp1(1)
Nel(2)   = N1(1)
k        = 2

!! Otherwise,  
!   Sp(1)  = sp1(1)
!   Nel(1) = N1(1)
!   k = 1

   num1(1) = k
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      if (Sp(k)==atomSp(4)) SpMass(k)=atomMass(4)
   print*, SpMass(k)

out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      if (Sp(k)==atomSp(4)) SpMass(k)=atomMass(4)
      print*, SpMass(k)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      if (Sp(k)==atomSp(4)) SpMass(k)=atomMass(4)
   end do out2
out3:do i=1,nb3
      do j=1,k
         if (sp3(i)==Sp(j)) then
            num3(i) = j
            cycle out3
         end if
      end do
      k = k+1
      num3(i) = k
      Sp(k) = sp3(i)
      Nel(k) = N3(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      if (Sp(k)==atomSp(4)) SpMass(k)=atomMass(4)
   end do out3


   !print*, a1, a2, a3
   !print*, a1Rot, a2Rot, a3Rot
!   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
!   write(u,'(a,i0)') 'NumberOfAtoms  ', (sz1+sz2+sz3)
!   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz2+sz3)
!   write(u,'(a)') 'LatticeConstant  1.00 Ang'
!   write(u,'(a)') '%block LatticeVectors'
!   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
!   write(u,'(a)') '%endblock LatticeVectors'
!   write(u,*)
!   write(u2,*)
   xlo = 0.0d0
   xhi = a1Rot(1)
   ylo = 0.0d0
   yhi = a2Rot(2)
   zlo = 0.0d0
   zhi = 35.0d0
   xy = a2Rot(1) - 0.000000001d0
   xz = 0.0d0
   yz = 0.0d0
   write(*,'(3f16.9)') (a1Rot(j),j=1,2)
   write(*,'(3f16.9)') (a2Rot(j),j=1,2)
!   write(u2,'(a,i0)') 'Lattice vectors:'
!   write(u2,'(a,3f16.9)') 'a:', xhi, 0.0d0, 0.0d0
!   write(u2,'(a,3f16.9)') 'b:', xy, yhi, 0.0d0
!   write(u2,'(a,3f16.9)') 'c:', 0.0d0, 0.0d0, 35.0d0
   !write(u2,'(a,i0)') 'Lattice vectors:'
   !do i=1,3
   !   write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   !end do
!   write(u2,*)
!   write(u2,'(a)') 'Fractional coordinates:'
!   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
!   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   write(u5,*)
   write(u5,'(a)') 'NNN atoms'
   write(u5,'(i1,a)') k, ' atom types'
   write(u5,*)
   write(u5,'(2f16.9,a)') xlo, xhi, ' xlo xhi'
   write(u5,'(2f16.9,a)') ylo, yhi, ' ylo yhi'
   write(u5,'(2f16.9,a)') zlo, zhi, ' zlo zhi'
   write(u5,'(3f16.9,a)') xy, xz, yz, ' xy xz yz'
   write(u5,*)
   write(u5,'(a)') ' Masses'
   write(u5,*)
   do j=1,k
   write(u5,'(i1, 3f16.9)') j, SpMass(j)
   write(*,'(i1, 3f16.9)') j, SpMass(j)
   end do   
   write(u5,*)
   write(u5,'(a)') ' Atoms'
   write(u5,*)
   j = 0
   rz = (z/2.0d0 - 3.0d0*h/2.0d0)/z
   write(*,'(f16.3)') rz*z
   !do i=1,sz1
   !   j = j+1
   !   !write(u,'(3f16.9,i8)') Xs1(:,i), rz, num1(mod(i+1,nb1)+1)
   !   !write(u2,'(A8,3f16.9)') sp1(mod(i+1,nb1)+1), Xs1(:,i), rz
   !   vecs(1)  = Rcell(1,1)*(Xs1(1,i))+Rcell(1,2)*(Xs1(2,i)) 
   !   vecs(2)  = Rcell(2,1)*(Xs1(1,i))+Rcell(2,2)*(Xs1(2,i)) 
   !   vecs(3)  = rz*z
   !   !write(u3,'(a,3f16.9,a)') 'basis ', Xs1(:,i), rz, ' &'
   !   !write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
   !   !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
   !   write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', num1(mod(i+1,nb1)+1), ' 0.0 ', vecs(:), 0,0,0 !Xs1(:,i), rz, 0, 0, 0 
   !   write(u6,'(i5,i8)') j,  mod(i+1,nb1)+1
   !   write(u7,'(i5,i8)') j,  1 
   !end do
   !rz1 = rz
   !   write(*,*)
   !   write(*,'(a)')       'Final N atoms  :'
   !   write(*,'(a,i8)')    '     L1        :', j
   j2 = 0
   rz = (z/2.0d0 - 1.0d0*h/2.0d0)/z
!   write(*,'(f16.3)') rz*z
   do i=1,sz2
      j = j +1
      j2= j2+1
      !write(u,'(3f16.9,i8)') Xs2(:,i), rz, num2(mod(i+1,nb2)+1)
      !write(u2,'(A8,3f16.9)') sp2(mod(i+1,nb2)+1), Xs2(:,i), rz
      vecs(1)  = Rcell(1,1)*(Xs2(1,i))+Rcell(1,2)*(Xs2(2,i)) 
      vecs(2)  = Rcell(2,1)*(Xs2(1,i))+Rcell(2,2)*(Xs2(2,i)) 
      vecs(3)  = rz*z
      !write(u3,'(a,3f16.9,a)') 'basis ', Xs2(:,i), rz, ' &'
      !write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', num2(mod(i+1,nb2)+1), ' 0.0 ', vecs(:), 0,0,0 ! Xs2(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j,  mod(i+1,nb2)+1+2
      write(u7,'(i5,i8)') j,  2 
   end do
   rz2 = rz
      write(*,'(a,i8)')    '     L2(ref)   :', j2
   j2 = 0
   j3 = 0
   rz = (z/2.0d0 + 1.0d0*h/2.0d0)/z
 do i=1,sz3
      f(1)  = ( Rcell0(1,1)*(Xs3(1,i))+Rcell0(1,2)*(Xs3(2,i)) -basisCirc0(1) )
      f(2)  = ( Rcell0(2,1)*(Xs3(1,i))+Rcell0(2,2)*(Xs3(2,i)) -basisCirc0(2) )

      if ((((    f(2)< sqrt(3.0d0)*(f(1)+a*(rradius+0.5d0)))  &
         &.and. (f(2)>-sqrt(3.0d0)*(f(1)+a*(rradius+0.5d0)))) &
         &.and.((f(2)<-sqrt(3.0d0)*(f(1)-a*(rradius+0.5d0)))  &
         &.and. (f(2)> sqrt(3.0d0)*(f(1)-a*(rradius+0.5d0)))))&
         &.and.((f(2)< sqrt(3.0d0)*0.5d0*a*(rradius+0.5d0))   &
         &.and. (f(2)>-sqrt(3.0d0)*0.5d0*a*(rradius+0.5d0)))) then

         j = j+1
         j2 = j2+1
         
         if (mod(i+1,nb3)+1 .eq. 1) then ! A sublattice
            f(1)  = delta32* ( Rcell(1,1)*(Xs3(1,i))+Rcell(1,2)*(Xs3(2,i))  &
                    & +ucell(1,1)*(basis32(1,1))+ucell(1,2)*(basis32(2,1)) -basisCirc(1) )
            f(2)  = delta32* ( Rcell(2,1)*(Xs3(1,i))+Rcell(2,2)*(Xs3(2,i))  &
                    & +ucell(2,1)*(basis32(1,1))+ucell(2,2)*(basis32(2,1)) -basisCirc(2) )
         else ! B sublattice
            f(1)  = delta32* ( Rcell(1,1)*(Xs3(1,i))+Rcell(1,2)*(Xs3(2,i))  &
                    & +ucell(1,1)*(basis32(1,2))+ucell(1,2)*(basis32(2,2)) -basisCirc(1) )
            f(2)  = delta32* ( Rcell(2,1)*(Xs3(1,i))+Rcell(2,2)*(Xs3(2,i))  &
                    & +ucell(2,1)*(basis32(1,2))+ucell(2,2)*(basis32(2,2)) -basisCirc(2) )
         end if

         vecs(1) = cos(phi32)*f(1) - sin(phi32)*f(2) + basisCirc(1)
         vecs(2) = sin(phi32)*f(1) + cos(phi32)*f(2) + basisCirc(2)
         vecs(3) = rz*z
         !write(u,'(3f16.9,i8)') Xs3(:,i)+bssCirc(:), rz, num3(mod(i+1,nb3)+1)
         !write(u,'(3f16.9,i8)') Xs3(:,i), rz, num3(mod(i+1,nb3)+1)
         !write(u2,'(A8,3f16.9)') sp3(mod(i+1,nb3)+1), Xs3(:,i), rz
         !write(u3,'(a,3f16.9,a)') 'basis ', Xs3(:,i), rz, ' &'
         !write(u4,'(a,i5,a,a)') 'basis ', j, ' 3', ' &'
         !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
         write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', num3(mod(i+1,nb3)+1), ' 0.0 ', vecs(:), 0,0,0 ! Xs3(:,i), rz, 0, 0, 0 
         !write(u6,'(i5,i8,i8)') j, 3, mod(i+1,nb3)+1
         write(u6,'(i5,i8)') j,  mod(i+1,nb3)+1+4
         write(u7,'(i5,i8)') j,  3 

      else  ! To add H passivation
         do nn = 1,3 ! For the three nearest A or B sublattices
         if (mod(i+1,nb3)+1 .eq. 1) then ! A sublattice
            if (nn.eq.1) then 
               nn1= 0.0d0; nn2= 0.0d0  
            else if (nn.eq.2) then 
               nn1=-1.0d0; nn2= 0.0d0 
            else if (nn.eq.3) then 
               nn1= 0.0d0; nn2=-1.0d0  
            end if
            f(1)  = ( Rcell0(1,1)*(Xs3(1,i+1))+Rcell0(1,2)*(Xs3(2,i+1)) -basisCirc0(1) ) + nn1*ucell0(1,1) + nn2*ucell0(1,2) 
            f(2)  = ( Rcell0(2,1)*(Xs3(1,i+1))+Rcell0(2,2)*(Xs3(2,i+1)) -basisCirc0(2) ) + nn1*ucell0(2,1) + nn2*ucell0(2,2) 
         else  ! B sublattice
            if (nn.eq.1) then 
               nn1= 0.0d0; nn2= 0.0d0  
            else if (nn.eq.2) then 
               nn1= 1.0d0; nn2= 0.0d0 
            else if (nn.eq.3) then 
               nn1= 0.0d0; nn2= 1.0d0  
            end if
            f(1)  = ( Rcell0(1,1)*(Xs3(1,i-1))+Rcell0(1,2)*(Xs3(2,i-1)) -basisCirc0(1) ) + nn1*ucell0(1,1) + nn2*ucell0(1,2) 
            f(2)  = ( Rcell0(2,1)*(Xs3(1,i-1))+Rcell0(2,2)*(Xs3(2,i-1)) -basisCirc0(2) ) + nn1*ucell0(2,1) + nn2*ucell0(2,2) 
         end if

         !dist2    = sqrt( vecs2(1)*vecs2(1) + vecs2(2)*vecs2(2) )

         if ((((    f(2)< sqrt(3.0d0)*(f(1)+a*(rradius+0.5d0)))  &
            &.and. (f(2)>-sqrt(3.0d0)*(f(1)+a*(rradius+0.5d0)))) &
            &.and.((f(2)<-sqrt(3.0d0)*(f(1)-a*(rradius+0.5d0)))  &
            &.and. (f(2)> sqrt(3.0d0)*(f(1)-a*(rradius+0.5d0)))))&
            &.and.((f(2)< sqrt(3.0d0)*0.5d0*a*(rradius+0.5d0))   &
            &.and. (f(2)>-sqrt(3.0d0)*0.5d0*a*(rradius+0.5d0)))) then
         !if (dist2 <= rradius*a) then
            j = j+1
            j3 = j3+1

            if (mod(i+1,nb3)+1 .eq. 1) then ! A sublattice
               f(1)  = delta32* ( Rcell(1,1)*(Xs3(1,i))+Rcell(1,2)*(Xs3(2,i))  &
                       & +ucell(1,1)*(basis32(1,1))+ucell(1,2)*(basis32(2,1)) -basisCirc(1) )
               f(2)  = delta32* ( Rcell(2,1)*(Xs3(1,i))+Rcell(2,2)*(Xs3(2,i))  &
                       & +ucell(2,1)*(basis32(1,1))+ucell(2,2)*(basis32(2,1)) -basisCirc(2) )
            else ! B sublattice
               f(1)  = delta32* ( Rcell(1,1)*(Xs3(1,i))+Rcell(1,2)*(Xs3(2,i))  &
                       & +ucell(1,1)*(basis32(1,2))+ucell(1,2)*(basis32(2,2)) -basisCirc(1) )
               f(2)  = delta32* ( Rcell(2,1)*(Xs3(1,i))+Rcell(2,2)*(Xs3(2,i))  &
                       & +ucell(2,1)*(basis32(1,2))+ucell(2,2)*(basis32(2,2)) -basisCirc(2) )
            end if

            vecs(1) = cos(phi32)*f(1) - sin(phi32)*f(2) + basisCirc(1)
            vecs(2) = sin(phi32)*f(1) + cos(phi32)*f(2) + basisCirc(2)
            vecs(3) = rz*z

            !write(u,'(3f16.9,i8)') Xs3(:,i), rz,1  !num1(mod(i+1,nb1)+1)
            !write(u2,'(A8,3f16.9)') '  H ' , Xs3(:,i), rz
            !write(u3,'(a,3f16.9,a)') 'basis ', Xs3(:,i), rz, ' &'
            !write(u4,'(a,i5,a,a)') 'basis ', j, ' 3', ' &'
            !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
            !write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 3 ', 1, ' 0.0 ', Xs3(:,i), rz, 0, 0, 0 
            write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', 1, ' 0.0 ', vecs(:), 0,0,0 ! Xs3(:,i), rz, 0, 0, 0 
            !write(u6,'(i5,i8,i8)') j, 3, mod(i+1,nb3)+1
            write(u6,'(i5,i8)') j,  mod(i+1,nb3)+1+6
            write(u7,'(i5,i8)') j,  4 
         end if
         end do
      end if
   end do
      write(*,'(a,i8)')    '     L3(Flake) :', j2
      write(*,'(a,i8)')    '     L4(H)     :', j3
      write(*,'(a,i8)')    '        Total  :', j
      write(*,*)
   write(*,'(f16.3)') rz1*z
   write(*,'(f16.3)') rz2*z
   write(*,'(f16.3)') rz*z
      write(*,*)
   !write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   !close(u)
   !close(u2)
   close(u5)
   close(u6)
   close(u7)

end subroutine MoireWrite3L_Hexa

subroutine MoireWrite2L_Circ(Xs1,Xs2,mm,a,h,z,fdfname,sp1,sp2,N1,N2,rradius,bssCirc,ucellCirc)

   implicit none

   integer, parameter :: u=10, u2=11, u3=12, u4=13, u5=14, u6=15

   double precision, intent(in) :: Xs1(:,:), Xs2(:,:),  a, h, z, rradius, bssCirc(:) , ucellCirc(:,:)
   integer, intent(in) :: mm(4), N1(:), N2(:)
   character(*), intent(in) :: fdfname, sp1(:), sp2(:)

   double precision :: ucell(2,2), ucellF(2,2) ,Rcell(3,3), m1(2,2), m2(2,2) !, m3(2,2) 
   double precision :: rz,rz1,ang
   integer :: sz1, sz2, i, j, k, j2, nb1, nb2, nn1,nn2,nn !, atomSp(3)=(/5,7,6/) !corresponding (/'B','N','C'/)
   character(60), pointer :: Sp(:)
   double precision, pointer :: SpMass(:)
   character :: atomSp(3)=(/'B','N','C'/)
   double precision :: atomMass(3)=(/10.800000191, 14.000000000, 12.010000229/)
   integer, pointer :: Nel(:),  num1(:), num2(:) 
   double precision :: phirad,phirad2, e1(3), dist, vecs(2), dist2, vecs2(2) , basisCirc(2)

   open(u,FILE=fdfname,STATUS='replace')
   open(u2,FILE="BLBL.frac",STATUS='replace')
   open(u3,FILE="BLBL.basis1",STATUS='replace')
   open(u4,FILE="BLBL.basis2",STATUS='replace')
   open(u5,FILE="BLBL.mol", STATUS='replace')
   open(u6,FILE="sublattices.dat", STATUS='replace')

   

   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)

   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
!   m3(:,1) = (/mm(5),mm(6)/)
!   m3(:,2) = (/-mm(6),mm(5)+mm(6)/)
   Rcell = 0.0d0
   Rcell(1:2,1:2) = matmul(ucell,m2)
   Rcell(3,3) = z

   ucellF = a*ucellCirc

   write(*,*)
   write(*,'(a,2f16.9)')   'Center in Frac.  : ', bssCirc(1), bssCirc(2)
   basisCirc(1) = Rcell(1,1)*bssCirc(1)+Rcell(1,2)*bssCirc(2)
   basisCirc(2) = Rcell(2,1)*bssCirc(1)+Rcell(2,2)*bssCirc(2)
   write(*,'(a,2f16.9)')   'Center in Cart.  : ', basisCirc(1), basisCirc(2)


      !f(2) = Rcell(1,1)*(X(2,i1)-X(1,i1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      !f(1) = (X(1,i1)-f(2)*Rcell(1,2))/Rcell(1,1)
      !X(:,i1) = mod(f,1.0d0)
      !if (X(1,i1)<0.0d0) X(1,i1) = X(1,i1) + 1.0d0
      !if (X(2,i1)<0.0d0) X(2,i1) = X(2,i1) + 1.0d0

!      ucellF(2,1) = Rcell(1,1)*(ucellCirc(2,1)-ucellCirc(1,1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
!      ucellF(1,1) = (ucellCirc(1,1)-ucellF(2,1)*Rcell(1,2))/Rcell(1,1)
!      ucellF(:,1) = mod(ucellF(:,1),1.0d0)
!      if (ucellF(1,1)<0.0d0) ucellF(1,1) = ucellF(1,1) + 1.0d0
!      if (ucellF(2,1)<0.0d0) ucellF(2,1) = ucellF(2,1) + 1.0d0
!      ucellF(2,2) = Rcell(1,1)*(ucellCirc(2,2)-ucellCirc(1,2)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
!      ucellF(1,2) = (ucellCirc(1,2)-ucellF(2,2)*Rcell(1,2))/Rcell(1,1)
!      ucellF(:,2) = mod(ucellF(:,2),1.0d0)
!      if (ucellF(1,2)<0.0d0) ucellF(1,2) = ucellF(1,2) + 1.0d0
!      if (ucellF(2,2)<0.0d0) ucellF(2,2) = ucellF(2,2) + 1.0d0

   sz1 = size(Xs1,2)
   sz2 = size(Xs2,2)
   nb1 = size(sp1)
   nb2 = size(sp2)

   allocate(Sp(nb1+nb2))
   allocate(Nel(nb1+nb2))
   allocate(SpMass(nb1+nb2))
   allocate(num1(nb1))
   allocate(num2(nb2))

   Sp(1)  = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = 1
   if (Sp(1)==atomSp(1)) SpMass(1)=atomMass(1)
   if (Sp(1)==atomSp(2)) SpMass(1)=atomMass(2)
   if (Sp(1)==atomSp(3)) SpMass(1)=atomMass(3)
   print*, SpMass(k)

out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      print*, SpMass(k)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
   end do out2

   write(u,'(a,i0)') 'NumberOfSpecies', k
   write(u,'(a)') '%block ChemicalSpeciesLabel'
   do i=1,k
      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
   end do
   !a1 = (/Rcell(1,1), Rcell(1,2), Rcell(1,3)/)
   !a2 = (/Rcell(2,1), Rcell(2,2), Rcell(2,3)/)
   !a3 = (/Rcell(3,1), Rcell(3,2), Rcell(3,3)/)
   a1 = Rcell(:,1)
   a2 = Rcell(:,2)
   a3 = Rcell(:,3)
   e1 = (/1,0,0/)
   phirad = acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
   write(*,'(f16.9, a)')  phirad*180.0d0/pi , ' : Angle of commensurate cell'
   Rot = reshape((/cos(-phirad),sin(-phirad),0.0d0,-sin(-phirad),cos(-phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
   a1Rot = matmul(Rot,a1)
   a2Rot = matmul(Rot,a2)
   a3Rot = matmul(Rot,a3)
   !print*, a1, a2, a3
   !print*, a1Rot, a2Rot, a3Rot
   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
   write(u,'(a,i0)') 'NumberOfAtoms  ', (sz1+sz2)
   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz2)
   write(u,'(a)') 'LatticeConstant  1.00 Ang'
   write(u,'(a)') '%block LatticeVectors'
   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
   write(u,'(a)') '%endblock LatticeVectors'
   write(u,*)
   write(u2,*)
   xlo = 0.0d0
   xhi = a1Rot(1)
   ylo = 0.0d0
   yhi = a2Rot(2)
   zlo = 0.0d0
   zhi = 35.0d0
   xy = a2Rot(1) - 0.000000001d0
   xz = 0.0d0
   yz = 0.0d0
   write(*,'(3f16.9)') (a1Rot(j),j=1,2)
   write(*,'(3f16.9)') (a2Rot(j),j=1,2)
   write(u2,'(a,i0)') 'Lattice vectors:'
   write(u2,'(a,3f16.9)') 'a:', xhi, 0.0d0, 0.0d0
   write(u2,'(a,3f16.9)') 'b:', xy, yhi, 0.0d0
   write(u2,'(a,3f16.9)') 'c:', 0.0d0, 0.0d0, 35.0d0
   !write(u2,'(a,i0)') 'Lattice vectors:'
   !do i=1,3
   !   write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   !end do
   write(u2,*)
   write(u2,'(a)') 'Fractional coordinates:'
   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   write(u5,*)
   write(u5,'(a)') 'NNN atoms'
   write(u5,'(i1,a)') k+1, ' atom types'
   write(u5,*)
   write(u5,'(2f16.9,a)') xlo, xhi, ' xlo xhi'
   write(u5,'(2f16.9,a)') ylo, yhi, ' ylo yhi'
   write(u5,'(2f16.9,a)') zlo, zhi, ' zlo zhi'
   write(u5,'(3f16.9,a)') xy, xz, yz, ' xy xz yz'
   write(u5,*)
   write(u5,'(a)') ' Masses'
   write(u5,*)
   do j=1,k
   write(u5,'(i1, 3f16.9)') j, SpMass(j)
   write(*,'(i1, 3f16.9)') j, SpMass(j)
   end do   
   write(u5,'(i1, 3f16.9)') k+1, 1.00784d0
   write(*,'(i1, 3f16.9)') k+1, 1.00784d0
   write(u5,*)
   write(u5,'(a)') ' Atoms'
   write(u5,*)
   j = 0
   rz = (z/2.0d0 + 1.0d0*h/2.0d0)/z
   write(*,'(f16.3)') rz*z
   do i=1,sz2
      j = j+1
      write(u,'(3f16.9,i8)') Xs2(:,i), rz, num2(mod(i+1,nb2)+1)
      write(u2,'(A8,3f16.9)') sp2(mod(i+1,nb2)+1), Xs2(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs2(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', num2(mod(i+1,nb2)+1), ' 0.0 ', Xs2(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8,i8)') j, 2, mod(i+1,nb2)+1+2
   end do
   rz1 = rz
      write(*,*)
      write(*,'(a)')       'Final N atoms  :'
      write(*,'(a,i8)')    '        Subs   :', j
   j2 = 0
   rz = (z/2.0d0 - 1.0d0*h/2.0d0)/z
 do i=1,sz1
      vecs(1)  = Rcell(1,1)*(Xs1(1,i))+Rcell(1,2)*(Xs1(2,i)) !-basisCirc(1)
      vecs(2)  = Rcell(2,1)*(Xs1(1,i))+Rcell(2,2)*(Xs1(2,i)) !-basisCirc(2)
      dist     = sqrt( vecs(1)*vecs(1) + vecs(2)*vecs(2) )
      if (dist <= rradius) then
         j = j+1
         j2 = j2+1
         write(u,'(3f16.9,i8)') Xs1(:,i)+bssCirc(:), rz, num1(mod(i+1,nb1)+1)
         write(u2,'(A8,3f16.9)') sp1(mod(i+1,nb1)+1), Xs1(:,i)+bssCirc(:), rz
         write(u3,'(a,3f16.9,a)') 'basis ', Xs1(:,i)+bssCirc(:), rz, ' &'
         write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
         !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
         write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', num1(mod(i+1,nb1)+1), ' 0.0 ', Xs1(:,i)+bssCirc(:), rz, 0, 0, 0 
         write(u6,'(i5,i8,i8)') j, 1, mod(i+1,nb1)+1

      else
         if (mod(i+1,nb1)+1 .eq. 1) then ! A sublattice
         do nn = 1,3 ! For the three nearest B sublattices
            if (nn.eq.1) then 
               nn1= 0.0d0; nn2= 0.0d0  
            else if (nn.eq.2) then 
               nn1=-1.0d0; nn2= 0.0d0 
            else if (nn.eq.3) then 
               nn1= 0.0d0; nn2=-1.0d0  
            end if
            vecs2(1)  = Rcell(1,1)*(Xs1(1,i+1))+Rcell(1,2)*(Xs1(2,i+1))+nn1*ucellF(1,1)+nn2*ucellF(1,2) !-basisCirc(1)
            vecs2(2)  = Rcell(2,1)*(Xs1(1,i+1))+Rcell(2,2)*(Xs1(2,i+1))+nn1*ucellF(2,1)+nn2*ucellF(2,2) !-basisCirc(2)
            dist2    = sqrt( vecs2(1)*vecs2(1) + vecs2(2)*vecs2(2) )
            if (dist2 <= rradius) then
               j = j+1
               j2 = j2+1
               write(u,'(3f16.9,i8)') Xs1(:,i)+bssCirc(:), rz,4  !num1(mod(i+1,nb1)+1)
               write(u2,'(A8,3f16.9)') '  H ' , Xs1(:,i)+bssCirc(:), rz
               write(u3,'(a,3f16.9,a)') 'basis ', Xs1(:,i)+bssCirc(:), rz, ' &'
               write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
               !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
               write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', 4, ' 0.0 ',  &
                                                 & Xs1(:,i)+bssCirc(:), rz, 0, 0, 0 
               write(u6,'(i5,i8,i8)') j, 1, mod(i+1,nb1)+1
            end if
         end do
         else ! B sublattice
         do nn = 1,3 ! For the three nearest A sublattices
            if (nn.eq.1) then 
               nn1= 0.0d0; nn2= 0.0d0  
            else if (nn.eq.2) then 
               nn1= 1.0d0; nn2= 0.0d0 
            else if (nn.eq.3) then 
               nn1= 0.0d0; nn2= 1.0d0  
            end if
            vecs2(1)  = Rcell(1,1)*(Xs1(1,i-1))+Rcell(1,2)*(Xs1(2,i-1))+nn1*ucellF(1,1)+nn2*ucellF(1,2) !-basisCirc(1)
            vecs2(2)  = Rcell(2,1)*(Xs1(1,i-1))+Rcell(2,2)*(Xs1(2,i-1))+nn1*ucellF(2,1)+nn2*ucellF(2,2) !-basisCirc(2)
            dist2    = sqrt( vecs2(1)*vecs2(1) + vecs2(2)*vecs2(2) )
            if (dist2 <= rradius) then
               j = j+1
               j2 = j2+1
               write(u,'(3f16.9,i8)') Xs1(:,i)+bssCirc(:), rz, 4!num1(mod(i+1,nb1)+1)
               write(u2,'(A8,3f16.9)') '  H ' , Xs1(:,i)+bssCirc(:), rz
               write(u3,'(a,3f16.9,a)') 'basis ', Xs1(:,i)+bssCirc(:), rz, ' &'
               write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
               !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
               write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', 4, ' 0.0 ', &
                                                  & Xs1(:,i)+bssCirc(:), rz, 0, 0, 0 
               write(u6,'(i5,i8,i8)') j, 1, mod(i+1,nb1)+1
            end if
         end do
         end if
      end if
   end do
      write(*,'(a,i8)')    '        Flake  :', j2
      write(*,'(a,i8)')    '        Total  :', j
      write(*,*)
   write(*,'(f16.3)') rz1*z
   write(*,'(f16.3)') rz*z
      write(*,*)
   write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   close(u)

end subroutine MoireWrite2L_Circ

subroutine MoireWrite3Leach(Xs1,Xs2,Xs3,mm,a,h,z,fdfname,sp1,sp2,sp3,N1,N2,N3)

   implicit none

   integer, parameter :: u=10, u2=11, u3=12, u4=13, u5=14, u6=15

   double precision, intent(in) :: Xs1(:,:), Xs2(:,:), Xs3(:,:),  a, h, z
   integer, intent(in) :: mm(6), N1(:), N2(:), N3(:)
   character(*), intent(in) :: fdfname, sp1(:), sp2(:), sp3(:)

   double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2), m3(2,2)
   double precision :: rz
   integer :: sz1, sz2, sz3, i, j, k, nb1, nb2, nb3 !, atomSp(3)=(/5,7,6/) !corresponding (/'B','N','C'/)
   character(60), pointer :: Sp(:)
   double precision, pointer :: SpMass(:)
   character :: atomSp(3)=(/'B','N','C'/)
   double precision :: atomMass(3)=(/10.800000191, 14.000000000, 12.010000229/)
   integer, pointer :: Nel(:),  num1(:), num2(:), num3(:)
   double precision :: phirad,phirad2, e1(3)

   open(u,FILE=fdfname,STATUS='replace')
   open(u2,FILE="BLBL.frac",STATUS='replace')
   open(u3,FILE="BLBL.basis1",STATUS='replace')
   open(u4,FILE="BLBL.basis2",STATUS='replace')
   open(u5,FILE="BLBL.mol", STATUS='replace')
   open(u6,FILE="sublattices.dat", STATUS='replace')

   

   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)

   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
   m3(:,1) = (/mm(5),mm(6)/)
   m3(:,2) = (/-mm(6),mm(5)+mm(6)/)
   Rcell = 0.0d0
   Rcell(1:2,1:2) = matmul(ucell,m2)
   Rcell(3,3) = z
   sz1 = size(Xs1,2)
   sz2 = size(Xs2,2)
   sz3 = size(Xs3,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
   nb3 = size(sp3)

   allocate(Sp(nb1+nb2+nb3))
   allocate(Nel(nb1+nb2+nb3))
   allocate(SpMass(nb1+nb2+nb3))
   allocate(num1(nb1))
   allocate(num2(nb2))
   allocate(num3(nb3))

   Sp(1)  = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = 1
   if (Sp(1)==atomSp(1)) SpMass(1)=atomMass(1)
   if (Sp(1)==atomSp(2)) SpMass(1)=atomMass(2)
   if (Sp(1)==atomSp(3)) SpMass(1)=atomMass(3)
   print*, SpMass(k)

out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
      print*, SpMass(k)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
   end do out2
out3:do i=1,nb3
      do j=1,k
         if (sp3(i)==Sp(j)) then
            num3(i) = j
            cycle out3
         end if
      end do
      k = k+1
      num3(i) = k
      Sp(k) = sp3(i)
      Nel(k) = N3(i)
      if (Sp(k)==atomSp(1)) SpMass(k)=atomMass(1)
      if (Sp(k)==atomSp(2)) SpMass(k)=atomMass(2)
      if (Sp(k)==atomSp(3)) SpMass(k)=atomMass(3)
   end do out3

   write(u,'(a,i0)') 'NumberOfSpecies', k
   write(u,'(a)') '%block ChemicalSpeciesLabel'
   do i=1,k
      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
   end do
   !a1 = (/Rcell(1,1), Rcell(1,2), Rcell(1,3)/)
   !a2 = (/Rcell(2,1), Rcell(2,2), Rcell(2,3)/)
   !a3 = (/Rcell(3,1), Rcell(3,2), Rcell(3,3)/)
   a1 = Rcell(:,1)
   a2 = Rcell(:,2)
   a3 = Rcell(:,3)
   e1 = (/1,0,0/)
   phirad = acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
   write(*,'(f16.9, a)')  phirad*180.0d0/pi , ' : Angle of commensurate cell'
   Rot = reshape((/cos(-phirad),sin(-phirad),0.0d0,-sin(-phirad),cos(-phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
   a1Rot = matmul(Rot,a1)
   a2Rot = matmul(Rot,a2)
   a3Rot = matmul(Rot,a3)
   !print*, a1, a2, a3
   !print*, a1Rot, a2Rot, a3Rot
   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
   write(u,'(a,i0)') 'NumberOfAtoms  ', (sz1+sz2+sz3)
   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz2+sz3)
   write(u,'(a)') 'LatticeConstant  1.00 Ang'
   write(u,'(a)') '%block LatticeVectors'
   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
   write(u,'(a)') '%endblock LatticeVectors'
   write(u,*)
   write(u2,*)
   xlo = 0.0d0
   xhi = a1Rot(1)
   ylo = 0.0d0
   yhi = a2Rot(2)
   zlo = 0.0d0
   zhi = 35.0d0
   xy = a2Rot(1) - 0.000000001d0
   xz = 0.0d0
   yz = 0.0d0
   write(*,'(3f16.9)') (a1Rot(j),j=1,2)
   write(*,'(3f16.9)') (a2Rot(j),j=1,2)
   write(u2,'(a,i0)') 'Lattice vectors:'
   write(u2,'(a,3f16.9)') 'a:', xhi, 0.0d0, 0.0d0
   write(u2,'(a,3f16.9)') 'b:', xy, yhi, 0.0d0
   write(u2,'(a,3f16.9)') 'c:', 0.0d0, 0.0d0, 35.0d0
   !write(u2,'(a,i0)') 'Lattice vectors:'
   !do i=1,3
   !   write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   !end do
   write(u2,*)
   write(u2,'(a)') 'Fractional coordinates:'
   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   write(u5,*)
   write(u5,'(i8,a)') (sz1+sz2+sz3), ' atoms'
   write(u5,'(i1,a)') k, ' atom types'
   write(u5,*)
   write(u5,'(2f16.9,a)') xlo, xhi, ' xlo xhi'
   write(u5,'(2f16.9,a)') ylo, yhi, ' ylo yhi'
   write(u5,'(2f16.9,a)') zlo, zhi, ' zlo zhi'
   write(u5,'(3f16.9,a)') xy, xz, yz, ' xy xz yz'
   write(u5,*)
   write(u5,'(a)') ' Masses'
   write(u5,*)
   do j=1,k
   write(u5,'(i1, 3f16.9)') j, SpMass(j)
   write(*,'(i1, 3f16.9)') j, SpMass(j)
   end do   
   write(u5,*)
   write(u5,'(a)') ' Atoms'
   write(u5,*)
   rz = (z/2.0d0 - 3.0d0*h/2.0d0)/z
   j = 0
   write(*,'(f16.3)') rz*z
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs1(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') sp1(mod(i+1,nb1)+1), Xs1(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs1(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', num1(mod(i+1,nb1)+1), ' 0.0 ', Xs1(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8,i8)') j, 1, mod(i+1,nb1)+1
   end do
   rz = (z/2.0d0 - 1.0d0*h/2.0d0)/z
   write(*,'(f16.3)') rz*z
   do i=1,sz2
      j = j+1
      write(u,'(3f16.9,i8)') Xs2(:,i), rz, num2(mod(i+1,nb2)+1)
      write(u2,'(A8,3f16.9)') sp2(mod(i+1,nb2)+1), Xs2(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs2(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', num2(mod(i+1,nb2)+1), ' 0.0 ', Xs2(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8,i8)') j, 2, mod(i+1,nb2)+1+2
   end do
   rz = (z/2.0d0 + 1.0d0*h/2.0d0)/z
   write(*,'(f16.3)') rz*z
   do i=1,sz3
      j = j+1
      write(u,'(3f16.9,i8)') Xs3(:,i), rz, num3(mod(i+1,nb3)+1)
      write(u2,'(A8,3f16.9)') sp3(mod(i+1,nb3)+1), Xs3(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs3(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', num3(mod(i+1,nb3)+1), ' 0.0 ', Xs3(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8,i8)') j, 3, mod(i+1,nb3)+1+4
   end do
   write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   close(u)

end subroutine MoireWrite3Leach

subroutine MoireWrite4L(Xs,Xs1,Xo,mm,a,h,z,name,sp1,sp2,N1,N2)

   implicit none

   integer, parameter :: u=10, u2=11, u3=12, u4=13, u5=14, u6=15

   double precision, intent(in) :: Xs(:,:), Xo(:,:), a, h, z
   double precision, intent(in) :: Xs1(:,:)
   integer, intent(in) :: mm(4), N1(:), N2(:)
   character(*), intent(in) :: name, sp1(:), sp2(:)

   double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2)
   double precision :: rz
   integer :: sz1, sz2, i, j, k, nb1, nb2
   character(60), pointer :: Sp(:)
   integer, pointer :: Nel(:), num1(:), num2(:)
   double precision :: phirad, e1(3)

   open(u,FILE=name,STATUS='replace')
   open(u2,FILE="BLBL.frac",STATUS='replace')
   open(u3,FILE="BLBL.basis1",STATUS='replace')
   open(u4,FILE="BLBL.basis2",STATUS='replace')
   open(u5,FILE="BLBL.mol", STATUS='replace')
   open(u6,FILE="sublattices.dat", STATUS='replace')
   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)
   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
   Rcell = 0.0d0
   Rcell(1:2,1:2) = matmul(ucell,m1)
   Rcell(3,3) = z
   sz1 = size(Xs,2)
   sz2 = size(Xo,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
   allocate(Sp(nb1+nb2))
   allocate(Nel(nb1+nb2))
   allocate(num1(nb1))
   allocate(num2(nb2))
   Sp(1) = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = 1
out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
   end do out2
   write(u,'(a,i0)') 'NumberOfSpecies', k
   write(u,'(a)') '%block ChemicalSpeciesLabel'
   do i=1,k
      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
   end do
   !a1 = (/Rcell(1,1), Rcell(1,2), Rcell(1,3)/)
   !a2 = (/Rcell(2,1), Rcell(2,2), Rcell(2,3)/)
   !a3 = (/Rcell(3,1), Rcell(3,2), Rcell(3,3)/)
   a1 = Rcell(:,1)
   a2 = Rcell(:,2)
   a3 = Rcell(:,3)
   e1 = (/1,0,0/)
   phirad = -acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
   Rot = reshape((/cos(phirad),sin(phirad),0.0d0,-sin(phirad),cos(phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
   a1Rot = matmul(Rot,a1)
   a2Rot = matmul(Rot,a2)
   a3Rot = matmul(Rot,a3)
   !print*, a1, a2, a3
   !print*, a1Rot, a2Rot, a3Rot
   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
   write(u,'(a,i0)') 'NumberOfAtoms  ', sz1+sz1+sz1+sz2
   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz1+sz1+sz2)
   write(u,'(a)') 'LatticeConstant  1.00 Ang'
   write(u,'(a)') '%block LatticeVectors'
   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
   write(u,'(a)') '%endblock LatticeVectors'
   write(u,*)
   write(u2,*)
   xlo = 0.0d0
   xhi = a1Rot(1)
   ylo = 0.0d0
   yhi = a2Rot(2)
   zlo = 0.0d0
   zhi = 35.0d0
   xy = a2Rot(1)
   xz = 0.0d0
   yz = 0.0d0
   write(u2,'(a,i0)') 'Lattice vectors:'
   write(u2,'(a,3f16.9)') 'a:', xhi, 0.0d0, 0.0d0
   write(u2,'(a,3f16.9)') 'b:', xy, yhi, 0.0d0
   write(u2,'(a,3f16.9)') 'c:', 0.0d0, 0.0d0, 35.0d0
   !write(u2,'(a,i0)') 'Lattice vectors:'
   !do i=1,3
   !   write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   !end do
   write(u2,*)
   write(u2,'(a)') 'Fractional coordinates:'
   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   write(u5,*)
   write(u5,'(i8,a)') sz1+sz1+sz1+sz2, ' atoms'
   write(u5,'(i1,a)') 2, ' atom types'
   write(u5,*)
   write(u5,'(2f16.9,a)') xlo, xhi, ' xlo xhi'
   write(u5,'(2f16.9,a)') ylo, yhi, ' ylo yhi'
   write(u5,'(2f16.9,a)') zlo, zhi, ' zlo zhi'
   write(u5,'(3f16.9,a)') xy, xz, yz, ' xy xz yz'
   write(u5,*)
   write(u5,'(a)') ' Masses'
   write(u5,*)
   write(u5,'(i1, 3f16.9)') 1, 12.0107
   write(u5,'(i1, 3f16.9)') 2, 12.0107
   write(u5,*)
   write(u5,'(a)') ' Atoms'
   write(u5,*)
   rz = (z/2.0d0 - h/2.0d0)/z
   j = 0
   print*, rz
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs1(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xs1(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs1(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', 1, ' 0.0 ', Xs1(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   rz = (z/2.0d0 + h/2.0d0)/z
   print*, rz
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs(:,i), rz, num2(mod(i+1,nb2)+1)
      write(u2,'(A8,3f16.9)') 'C', Xs(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', 1, ' 0.0 ', Xs(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   rz = (z/2.0d0 + 3.0*h/2.0d0)/z
   print*, rz
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs1(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xs1(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs1(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 2 ', 1, ' 0.0 ', Xs1(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   rz = (z/2.0d0 + 5.0*h/2.0d0)/z
   print*, rz
   do i=1,sz2
      j = j+1
      write(u,'(3f16.9,i8)') Xo(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xo(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xo(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      if (num1(mod(i+1,nb2)+1).eq.1) then
         write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', 2, ' 0.0 ', Xo(:,i), rz, 0, 0, 0 
      else
         write(u5,'(i8,a,i6,a,3f16.9,3i3)') j, ' 1 ', 3, ' 0.0 ', Xo(:,i), rz, 0, 0, 0 
      end if
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   close(u)

end subroutine MoireWrite4L

function gcd(i1,i2)

    integer :: gcd
    integer, intent(in) :: i1,i2
    integer :: a,b

    a = i1
    b = i2
    if (a > b) then
        gcd = a
        a = b
        b = gcd
    end if

    do
      gcd = mod(a, b)
      if (gcd == 0) exit
      a = b
      b = gcd
    end do

    gcd = b

end function gcd

end program GenMoire
