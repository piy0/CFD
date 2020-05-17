!=======================================================================================================================
!==== [ Module Definition ] ============================================================================================
!=======================================================================================================================
Module Definition
   Implicit none

   !--- Input Parameter ---
   Real*8 ,Parameter :: Lx = 10.0d0
   Real*8 ,Parameter :: Ly = 1.0d0
   Integer,Parameter :: Nx = 151
   Integer,Parameter :: Ny = 51
   Integer,Parameter :: ltime = 1000

   Real*8 ,Parameter :: dt = 0.01d0
   Real*8 ,Parameter :: Re = 200.0d0
   Real*8 ,Parameter :: r1_dt = 1.0d0/dt
   Real*8 ,Parameter :: r1_Re = 1.0d0/Re
   !-----------------------

   !--- Coordinate ---
   Real*8 x(0:Nx)
   Real*8 y(0:Ny)
   Real*8 dx,dy

   !--- variable ---
   Real*8 U(0:Nx,0:Ny)
   Real*8 V(0:Nx,0:Ny)
   Real*8 P(0:Nx,0:Ny)

   !--- differential ---
   Real*8 dU(1:Nx-1,1:Ny-1)
   Real*8 dV(1:Nx-1,1:Ny-1)
End Module

!=======================================================================================================================
!==== [ Program Main ] =================================================================================================
!=======================================================================================================================
Program Main
   use Definition
   Implicit none


   Call channel_flow


   Stop
End Program

!=======================================================================================================================
!==== [ Subroutine channel_flow ] ======================================================================================
!=======================================================================================================================
Subroutine channel_flow
   use Definition
   Implicit none

   Integer i,j,k,t

   Call Grid_Gene
   Call Ini_Condition
   Call Bou_Condition

   Do t=1,ltime
      Write(*,*) t

      Call Calc_NSE_term
      Call Calc_Velocity
      Call Calc_Poisson
   End Do

   open(100,file='100_output.txt')
   open(101,file='101_output.txt')
   Do j=0,Ny
   Do i=0,Nx
      Write(100,*) x(i),y(j),U(i,j)
      If(i.EQ.Nx-1) Write(101,*) y(j),U(i,j)
   End Do ; End Do
   close(100)
   close(101)

   Return
End Subroutine

!=======================================================================================================================
!==== [ Subroutine Grid_Gene ] =========================================================================================
!=======================================================================================================================
Subroutine Grid_Gene
   use Definition
   Implicit none

   Integer i,j,k

   !x-coordinate
   x = 0
   dx = Lx / Real(Nx)
   Do i=1,Nx
      x(i) = x(i-1) + dx
   End Do

   !y-coordinate
   y = 0
   dy = Ly / Real(Ny)
   Do j=1,Ny
      y(j) = y(j-1) + dy
   End Do

   Return
End Subroutine

!=======================================================================================================================
!==== [ Subroutine Ini_Condition ] =====================================================================================
!=======================================================================================================================
Subroutine Ini_Condition
   use Definition
   Implicit none

   Integer i,j,k

   !--- Initial_Condition ---
   Do j=0,Ny
   Do i=0,Nx
      U(i,j) = 1.0d0
      V(i,j) = 0.0d0
      P(i,j) = 0.0d0
   End Do ; End Do

   !--- output ---
   open(14,file='14_Ini_Condition.txt')
   Do j=0,Ny
   Do i=0,Nx
      Write(14,*) x(i),y(j),U(i,j),V(i,j),P(i,j)
   End Do ; End Do
   close(14)
   !---------------

   Return
End Subroutine

!=======================================================================================================================
!==== [ Subroutine Bou_Condition ] =====================================================================================
!=======================================================================================================================
Subroutine Bou_Condition
   use Definition
   Implicit none

   Integer i,j,k

   !--- InFlow ---
   Do j=0,Ny
      U(0,j) = 1.0d0
      V(0,j) = 0.0d0
      P(0,j) = P(1,j)
   End Do

   !--- OutFlow ---
   Do j=0,Ny
      U(Nx,j) = U(Nx-1,j)
      V(Nx,j) = V(Nx-1,j)
      P(Nx,j) = P(Nx-1,j)
   End Do

   !--- Wall ---
   Do i=0,Nx
      U(i,0) = 0.0d0
      V(i,0) = 0.0d0
      P(i,0) = P(i,1)

      U(i,Ny) = 0.0d0
      V(i,Ny) = 0.0d0
      P(i,Ny) = P(i,Ny-1)
   End Do

   !--- output ---
   open(15,file='15_Bou_Condition.txt')
   Do j=0,Ny
   Do i=0,Nx
      Write(15,*) x(i),y(j),U(i,j),V(i,j)
   End Do ; End Do
   close(15)
   !---------------

   Return
End Subroutine

!=======================================================================================================================
!==== [ Subroutine Calc_NSE_term ] =====================================================================================
!=======================================================================================================================
Subroutine Calc_NSE_term
   use Definition
   Implicit none

   Integer i,j,k

   Real*8 Ux(1:Nx-1,1:Ny-1),Uy(1:Nx-1,1:Ny-1),Uxx(1:Nx-1,1:Ny-1),Uyy(1:Nx-1,1:Ny-1)
   Real*8 Vx(1:Nx-1,1:Ny-1),Vy(1:Nx-1,1:Ny-1),Vxx(1:Nx-1,1:Ny-1),Vyy(1:Nx-1,1:Ny-1)

   Do j=1,Ny-1
   Do i=1,Nx-1
      Ux (i,j) = ( U(i+1,j )-U(i-1,j ) ) / (2.0d0*dx)
      Uxx(i,j) = ( U(i+1,j )-2.0d0*U(i,j)+U(i-1,j ) )/(dx*dx)

      Uy (i,j) = ( U(i ,j+1)-U(i ,j-1) ) / (2.0d0*dy)
      Vx (i,j) = ( V(i+1,j )-V(i-1,j ) ) / (2.0d0*dx)
      Vy (i,j) = ( V(i ,j+1)-V(i ,j-1) ) / (2.0d0*dy)
      Uyy(i,j) = ( U(i ,j+1)-2.0d0*U(i,j)+U(i ,j-1) )/(dy*dy)
      Vxx(i,j) = ( V(i+1,j )-2.0d0*V(i,j)+V(i-1,j ) )/(dx*dx)
      Vyy(i,j) = ( V(i ,j+1)-2.0d0*V(i,j)+V(i ,j-1) )/(dy*dy)

      dU(i,j) = r1_Re*( Uxx(i,j)+Uyy(i,j) ) - ( U(i,j)*Ux(i,j)+V(i,j)*Uy(i,j) )
      dV(i,j) = r1_Re*( Vxx(i,j)+Vyy(i,j) ) - ( U(i,j)*Vx(i,j)+V(i,j)*Vy(i,j) )
   End Do ; End Do

   Return
End Subroutine

!=======================================================================================================================
!==== [ Subroutine Calc_Velocity ] =====================================================================================
!=======================================================================================================================
Subroutine Calc_Velocity
   use Definition
   Implicit none

   Integer i,j,k

   Do j=1,Ny-1
   Do i=1,Nx-1
      U(i,j) = U(i,j) + dU(i,j)*dt
      V(i,j) = V(i,j) + dV(i,j)*dt
   End Do ; End Do

   Call Bou_Condition

   Return
End Subroutine

!=======================================================================================================================
!==== [ Subroutine Calc_Poisson ] ======================================================================================
!=======================================================================================================================
Subroutine Calc_Poisson
   use Definition
   Implicit none

   Integer i,j,k,l

   Real*8 P_prev,MaxErr,CurErr,f(0:Nx,0:Ny)
   Real*8 :: MaxP = 1.0d-10
   Real*8 :: Conv = 1.0d-3
   Integer,Parameter :: lmax=5000

   !--- Calc Poisson-Equation ---
   Do l=1,lmax

      MaxErr = 0.0d0
      CurErr = 0.0d0
      Do j=1,Ny-1
      Do i=1,Nx-1
         f(i,j) = r1_dt*(   ( U(i+1,j  )-U(i-1,j  ) ) / (2.0d0*dx)   &
                          + ( V(i  ,j+1)-V(i  ,j-1) ) / (2.0d0*dy) )
         P_prev = P(i,j)
         P(i,j) = (  dx*dx*(P(i,j+1)+P(i,j-1))                       &
                    +dy*dy*(P(i+1,j)+P(i-1,j))                       &
                    -dx*dx*dy*dy*f(i,j) ) / (2.0d0*(dx*dx+dy*dy))
         If (MaxP.LT.abs(P(i,j))) MaxP = P(i,j)
         CurErr = ( abs(P(i,j)-P_prev) ) / MaxP

         If (MaxErr.LT.CurErr) MaxErr = CurErr
      End Do ; End Do

      If (MaxErr.LT.Conv) Exit
   End Do

   Do j=1,Ny-1
   Do i=1,Nx-1
      U(i,j) = U(i,j) - dt*( P(i+1,j )-P(i-1,j ) ) / (2.0d0*dx)
      V(i,j) = V(i,j) - dt*( P(i ,j+1)-P(i ,j-1) ) / (2.0d0*dy)
   End Do ; End Do

   Call Bou_Condition

   Return
End Subroutine
