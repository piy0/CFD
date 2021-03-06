!=======================================================================================================================
!==== [ Module Definition ] ============================================================================================
!=======================================================================================================================
Module Definition

   Implicit none

   Real*8 ,Parameter :: pi = 3.14159265d0

   !=== Input Parameter ===
   Real*8 ,Parameter :: Lx = 25.0d0
   Real*8 ,Parameter :: Ly = 15.0d0
   Integer,Parameter :: Nx = 251
   Integer,Parameter :: Ny = 151
   Integer,Parameter :: ltime = 7000

   Real*8 ,Parameter :: dt = 0.01d0
   Real*8 ,Parameter :: Re = 150.0d0
   Real*8 ,Parameter :: r1_dt = 1.0d0/dt
   Real*8 ,Parameter :: r1_Re = 1.0d0/Re

   !--- Calc-Range ---
   Real*8 ,Parameter :: Lx_sta = Lx/5.0d0 - 0.5d0
   Real*8 ,Parameter :: Lx_end = Lx/5.0d0 + 0.5d0
   Real*8 ,Parameter :: Ly_sta = Ly/2.0d0 - 0.5d0
   Real*8 ,Parameter :: Ly_end = Ly/2.0d0 + 0.5d0
   Integer i_Lxsta, i_Lxend
   Integer j_Lysta, j_Lyend
   !------------------
   !=======================

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
Program main
   use Definition
   Implicit none

   Call backstep_flow

   Stop
End Program


!=======================================================================================================================
!==== [ Subroutine backstep_flow ] =====================================================================================
!=======================================================================================================================
Subroutine backstep_flow
   use Definition
   Implicit none

   Integer t

   Call Grid_Gene
   Call Calc_range
   Call Ini_Condition
   Call Bou_Condition

   Call output_vtk(0)

   !stop 

   Do t=1,ltime
      Write(*,*) t

      Call Calc_NSE_term
      Call Calc_Velocity
      Call Calc_Poisson

      If ( mod(t,10)==0 ) Call output_vtk(t)
   End Do

   Call output

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

   !--- output ----
   open(10,file='10_xGrid.txt')
   Do j=0,Ny
   Do i=0,Nx
      Write(10,*) x(i),y(j)
      If (i.EQ.Nx) Write(10,*)
   End Do ; End Do
   close(10)

   open(11,file='11_yGrid.txt')
   Do i=0,Nx
   Do j=0,Ny
      Write(11,*) x(i),y(j)
      If (j.EQ.Ny) Write(11,*)
   End Do ; End Do
   close(11)

   open(12,file='12_xyGrid.txt')
   Do j=0,Ny
   Do i=0,Nx
      Write(12,*) x(i),y(j)
   End Do ; End Do
   close(12)
   !---------------

   Return
End Subroutine


!=======================================================================================================================
!==== [ Subroutine Calc_range ] ========================================================================================
!=======================================================================================================================
Subroutine Calc_range
   use Definition
   Implicit none

   Integer i,j,k

   !x-range
   Do i=0,Nx
   	  If ( Lx_sta <= x(i) ) Then
         i_Lxsta = i
         Exit
      End If
   End Do

   Do i=0,Nx
   	  If ( Lx_end <= x(i) ) Then
         i_Lxend = i
         Exit
      End If
   End Do


   !y-range
   Do j=0,Ny
   	  If ( Ly_sta <= y(j) ) Then
   	  	 j_Lysta = j
   	     Exit
   	  End If
   End Do

   Do j=0,Ny
   	  If ( Ly_end <= y(j) ) Then
   	  	 j_Lyend = j
   	     Exit
   	  End If
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


   !--- Left ---
   Do j=0,Ny
      U(0,j) = 1.0d0
      V(0,j) = 0.0d0
      P(0,j) = P(1,j)
   End Do

   !--- Right ---
   Do j=0,Ny
      U(Nx,j) = U(Nx-1,j)
      V(Nx,j) = V(Nx-1,j)
      P(Nx,j) = P(Nx-1,j)
   End Do

   !--- Up ---
   Do i=0,Nx
      U(i,Ny) = 1.0d0
      V(i,Ny) = 0.0d0
      P(i,Ny) = P(i,Ny-1)
   End Do

   !--- Down ---
   Do i=0,Nx
      U(i,0) = 1.0d0
      V(i,0) = 0.0d0
      P(i,0) = P(i,1)
   End Do


   !--- in object ---
   Do j=j_Lysta,j_Lyend ; Do i=i_Lxsta,i_Lxend
      U(i,j) = 0.0d0
      V(i,j) = 0.0d0
   End Do ; End Do

   Do j=j_Lysta,j_Lyend
   	  P(i_Lxsta,j) = P(i_Lxsta+1,j)
   End Do

   Do i=i_Lxsta,i_Lxend
   	  P(i,j_Lysta) = P(i,j_Lysta+1)
   End Do


   !!--- output ---
   !open(15,file='15_Bou_Condition.txt')
   !Do j=0,Ny
   !Do i=0,Nx
   !   Write(15,*) x(i),y(j),U(i,j),V(i,j)
   !End Do ; End Do
   !close(15)
   !!---------------

   Return
End Subroutine


!=======================================================================================================================
!==== [ Subroutine Calc_NSE_term ] =====================================================================================
!=======================================================================================================================
Subroutine Calc_NSE_term
   use Definition
   Implicit none

   Integer i,j,k

   Real*8  Ux(1:Nx-1,1:Ny-1),Uy(1:Nx-1,1:Ny-1),Uxx(1:Nx-1,1:Ny-1),Uyy(1:Nx-1,1:Ny-1)
   Real*8  Vx(1:Nx-1,1:Ny-1),Vy(1:Nx-1,1:Ny-1),Vxx(1:Nx-1,1:Ny-1),Vyy(1:Nx-1,1:Ny-1)


   Do j=1,Ny-1
   Do i=1,Nx-1
      Ux (i,j) = ( U(i+1,j  )-U(i-1,j  ) ) / (2.0d0*dx)
      Uxx(i,j) = ( U(i+1,j  )-2.0d0*U(i,j)+U(i-1,j  ) )/(dx*dx)


      Uy (i,j) = ( U(i  ,j+1)-U(i  ,j-1) ) / (2.0d0*dy)
      Vx (i,j) = ( V(i+1,j  )-V(i-1,j  ) ) / (2.0d0*dx)
      Vy (i,j) = ( V(i  ,j+1)-V(i  ,j-1) ) / (2.0d0*dy)
      Uyy(i,j) = ( U(i  ,j+1)-2.0d0*U(i,j)+U(i  ,j-1) )/(dy*dy)
      Vxx(i,j) = ( V(i+1,j  )-2.0d0*V(i,j)+V(i-1,j  ) )/(dx*dx)
      Vyy(i,j) = ( V(i  ,j+1)-2.0d0*V(i,j)+V(i  ,j-1) )/(dy*dy)

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
   Real*8 :: Conv = 1.0d-5
   Integer,Parameter :: lmax=10000


   
   !--- Calc Velocity ---
   Do l=1,lmax
      
      MaxErr = 0.0d0
      CurErr = 0.0d0
      Do j=1,Ny-1
      Do i=1,Nx-1
         f(i,j) = r1_dt*( ( U(i+1,j  )-U(i-1,j  ) ) / (2.0d0*dx) + ( V(i  ,j+1)-V(i  ,j-1) ) / (2.0d0*dy) )
         P_prev = P(i,j)
         P(i,j) = ( dx*dx*(P(i,j+1)+P(i,j-1))+dy*dy*(P(i+1,j)+P(i-1,j))-dx*dx*dy*dy*f(i,j) ) / (2.0d0*(dx*dx+dy*dy))

         If (MaxP.LT.abs(P(i,j))) MaxP = P(i,j)
         CurErr = ( abs(P(i,j)-P_prev) ) / MaxP

         If (MaxErr.LT.CurErr) MaxErr = CurErr
      End Do ; End Do

      If (MaxErr.LT.Conv) Exit
   End Do

   Do j=1,Ny-1
   Do i=1,Nx-1
      U(i,j) = U(i,j) - dt*( P(i+1,j  )-P(i-1,j  ) ) / (2.0d0*dx)
      V(i,j) = V(i,j) - dt*( P(i  ,j+1)-P(i  ,j-1) ) / (2.0d0*dy)
   End Do ; End Do

   Call Bou_Condition

   Return
End Subroutine


!=======================================================================================================================
!==== [ Subroutine output ] ============================================================================================
!=======================================================================================================================
Subroutine output
   use Definition
   Implicit none

   Integer i,j,k


   open(100,file='100_output.txt')
   open(101,file='101_output.txt')
   Do j=0,Ny
   Do i=0,Nx
      Write(100,*) x(i),y(j),U(i,j)
      If(i.EQ.(Nx-1)/2) Write(101,*) y(j),U(i,j)
   End Do ; End Do
   close(100)
   close(101)


   Return
End Subroutine


!=======================================================================================================================
!==== [ Subroutine output_vtk ] ========================================================================================
!=======================================================================================================================
Subroutine output_vtk(step)
   use Definition
   Implicit none

   Integer i,j,k

   Integer step, points
   character(len=40) filename
  
   Write(filename,"(a,i5.5,a)") "output/output",int(step),".vtk"
   Open(10,file=filename)
   Write(10,"('# vtk DataFile Version 3.0')")
   Write(10,"('test')")
   Write(10,"('ASCII ')")
 
   Write(10,"('DATASET STRUCTURED_GRID')")
   Write(10,"('DIMENSIONS ',3(1x,i3))") Nx+1, Ny+1, 1
 
   Write(10,"('POINTS ',i9,' float')") (Nx+1)*(Ny+1)*1
   Do j=0,Ny ; Do i=0,Nx
      Write(10,"(3(f9.4,1x))") i*dx, j*dy, 0.0d0
   End Do ; End Do
 
   Write(10,"('POINT_DATA ',i9)") (Nx+1)*(Ny+1)*1
 
   !date input
   Write(10,"('SCALARS U float')")
   Write(10,"('LOOKUP_TABLE default')")
   Do j=0,Ny ; Do i=0,Nx
   	  Write(10,*)  U(i,j)
   End Do ; End Do


   Close(10)

   Return
End Subroutine