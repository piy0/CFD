!=======================================================================================================================
!========= [ Module Definition ] =======================================================================================
!=======================================================================================================================
Module Definition


   Implicit none

   Real*8,Parameter :: pi  = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825d0

   !=== Parameters ====================================================
   Integer,Parameter :: Nx = 201
   Integer,Parameter :: Ny = 201
   Integer,Parameter :: Nz = 101
   Real*8,Parameter  :: Lx = 10.0d0
   Real*8,Parameter  :: Ly = 10.0d0
   Real*8,Parameter  :: Lz = 10.0d0
   Real*8,Parameter  :: Re = 1000.0d0
   Real*8,Parameter  :: Kd = 5.3d0
   Real*8,Parameter  :: Fr = 0.6d0
   Real*8,Parameter  :: Pr = 1.0d0

   Integer,Parameter :: nlast = 100

   Real*8,Parameter  :: alpha_x  = 1.5d0
   Real*8,Parameter  :: alpha_y  = 2.5d0
   Real*8,Parameter  :: eps = 0.01d0
   !-------------------------------------------------------------------


   !=== Set_Parameter ===
   Real*8 x(-2:Nx+2), y(-2:Ny+2), z(-2:Nz+2), dz
   Real*8 dx(-1:Nx+2), dxp(-1:Nx+2), dxm(-1:Nx+2), c1_dx(-1:Nx+2), c1_dx2(-1:Nx+2), c1_dxm(-1:Nx+2), c1_dxp(-1:Nx+2)
   Real*8 dy(-1:Ny+2), dyp(-1:Ny+2), dym(-1:Ny+2), c1_dy(-1:Ny+2), c1_dy2(-1:Ny+2), c1_dym(-1:Ny+2), c1_dyp(-1:Ny+2)
   Real*8 dxhm(0:Nx), dxhp(0:Nx), c1_dxhm(0:Nx), c1_dxhp(0:Nx)
   Real*8 dyhm(0:Ny), dyhp(0:Ny), c1_dyhm(0:Ny), c1_dyhp(0:Ny)
   Real*8 c1_dxRe (0:Nx), c1_dxRePr(0:Nx), c1_dxhpRe(0:Nx)
   Real*8 c1_dyRe (0:Ny), c1_dyRePr(0:Ny), c1_dyhpRe(0:Ny)
   Real*8 c1_dxhmx(0:Nx), c1_dxhpx(0:Nx)
   Real*8 c1_dyhmy(0:Ny), c1_dyhpy(0:Ny)

   Real*8 c1_dtime, c1_pi, c1_Re, c1_RePr
   Real*8 c1_Nx, c1_Ny , c1_Nz  , c1_NxNy, c1_NxNz
   Real*8 c1_dz, c1_dz2, c1_dzRe, c1_dz2Re
   Real*8 c1_dzRePr,  c1_dz2RePr, c1_Fr2
   Real*8,Parameter :: c1_Fr =  1.0d0 / Fr


   !=== Set_IBM ===
   Real*4  x_U(-2:Nx+2),x_V(-2:Nx+2),x_W(-2:Nx+2),x_Ca(-2:Nx+2)
   Real*4  y_U(-2:Ny+2),y_V(-2:Ny+2),y_W(-2:Ny+2),y_Ca(-2:Ny+2)
   Real*4  z_U(-2:Nz+2),z_V(-2:Nz+2),z_W(-2:Nz+2),z_Ca(-2:Nz+2)

   !=== Generate_Object ===
   Integer i_OL,i_OR
   Integer j_OU,j_OD
   Real*4 ,Parameter :: x_Ob0 = Lx / 2.0d0 !中心x座標
   Real*4 ,Parameter :: y_Ob0 = Ly / 2.0d0 !中心y座標
   Real*4 ,Parameter :: Lx_Ob = 0.5d0
   Real*4 ,Parameter :: Ly_Ob = 0.5d0
   Real*4 ,Parameter :: ch_L  = 1.0d0 / Lx_Ob

   Real*4  x_ObU(-2:Nx+2),x_ObV (-2:Nx+2),y_ObUV (-2:Nx+2),y_ObUU(-2:Nx+2),y_ObDV (-2:Nx+2),y_ObDU(-2:Nx+2)
   Real*4  x_ObW(-2:Nx+2),x_ObCa(-2:Nx+2),y_ObUCa(-2:Nx+2),y_ObUW(-2:Nx+2),y_ObDCa(-2:Nx+2),y_ObDW(-2:Nx+2)

   !=== Mark_BG ===
   Integer j_ObUU(-2:Nx+2),j_ObUV(-2:Nx+2),j_ObUW(-2:Nx+2),j_ObUCa(-2:Nx+2)
   Integer j_ObDU(-2:Nx+2),j_ObDV(-2:Nx+2),j_ObDW(-2:Nx+2),j_ObDCa(-2:Nx+2)
   Integer IB_U(-2:Nx+2,-2:Ny+2),IB_V(-2:Nx+2,-2:Ny+2),IB_W(-2:Nx+2,-2:Ny+2),IB_Ca(-2:Nx+2,-2:Ny+2)
   Integer FS_U(-2:Nx+2,-2:Ny+2),FS_V(-2:Nx+2,-2:Ny+2),FS_W(-2:Nx+2,-2:Ny+2),FS_Ca(-2:Nx+2,-2:Ny+2)

   !=== Mark_BG ===
   Real*4 x_ObOU(-2:Nx+2),y_ObOUU(-2:Nx+2),y_ObODU(-2:Nx+2),x_IPU(-2:Nx+2),y_IPUU(-2:Nx+2),y_IPDU(-2:Nx+2)



   Real*8,Parameter :: c1_2      = 1.0d0 /   2.0d0
   Real*8,Parameter :: c1_3      = 1.0d0 /   3.0d0
   Real*8,Parameter :: c2_3      = 2.0d0 /   3.0d0
   Real*8,Parameter :: c3_2      = 3.0d0 /   2.0d0
   Real*8,Parameter :: c1_4      = 1.0d0 /   4.0d0
   Real*8,Parameter :: c1_5      = 1.0d0 /   5.0d0
   Real*8,Parameter :: c1_6      = 1.0d0 /   6.0d0
   Real*8,Parameter :: c2_6      = 2.0d0 /   6.0d0
   Real*8,Parameter :: c4_3      = 4.0d0 /   3.0d0
   Real*8,Parameter :: c1_8      = 1.0d0 /   8.0d0
   Real*8,Parameter :: c5_8      = 5.0d0 /   8.0d0
   Real*8,Parameter :: c9_8      = 9.0d0 /   8.0d0
   Real*8,Parameter :: c1_12     = 1.0d0 /  12.0d0
   Real*8,Parameter :: c1_100    = 1.0d0 / 100.0d0


End Module



!=======================================================================================================================
!========= [ Program Main ] ============================================================================================
!=======================================================================================================================
Program Main

   Use Definition
   Implicit none 
   Integer i, j, k



   Call Calc_CBW



   Stop
End Program


!=======================================================================================================================
!========= [ Subroutine Calc_CBW ] =====================================================================================
!=======================================================================================================================
Subroutine Calc_CBW

   Use Definition
   Implicit none 
   Integer i, j, k ,t


   !--- Set Coordinates and Parameters ---
   Call Set_Parameter
   Call Set_IBM


   Return
End Subroutine


!=======================================================================================================================
!========= [ Subroutine Set_Parameter ] ================================================================================
!=======================================================================================================================
Subroutine Set_Parameter

   Use Definition
   Implicit none 
   Integer i, j, k


   x = 0.0d0
   y = 0.0d0
   z = 0.0d0

   !--- x ---
   Do i = 0, Nx+2, 1
     !x(i) = x(i-1) + Lx/Real(Nx+2)!
     x(i) = 0.5d0*Lx*(1.0d0-atanh(tanh(alpha_x)*(1.0d0-2.0d0*(Dble(i)/Dble(Nx-1)))) / alpha_x )!x(i-1) + Lx/Real(Nx+2)
   End Do

   Do i = 1, Nx+2, 1
     dx(i) = x(i) - x(i-1)  
   End Do
     dx(0) = x(0) - ( 2.0d0 * x(0) - x(1) )

   Do i = 0, Nx, 1
     dxm(i)    = 0.5d0 * dx(i)
     dxp(i)    = 0.5d0 * dx(i)

     c1_dx (i) = 1.0d0 /  dx(i)
     c1_dx2(i) = 1.0d0 / (dx(i  ) * dx(i))

     c1_dxm(i) = 1.0d0 / dxm (i)
     c1_dxp(i) = 1.0d0 / dxp (i)
   End Do

   Do i = 1, Nx-1, 1
     dxhm(i)     = ( dx(i-1) + dx(i  ) ) * 0.5d0
     dxhp(i)     = ( dx(i  ) + dx(i+1) ) * 0.5d0

     c1_dxhm(i)  = 1.0d0 /   dxhm(i)
     c1_dxhp(i)  = 1.0d0 /   dxhp(i)
   End Do

   Do i = 1, Nx-1, 1
     c1_dxRe(i)   = c1_dx  (i) * c1_Re
     c1_dxhpRe(i) = c1_dxhp(i) * c1_Re
     c1_dxRePr(i) = c1_dx  (i) * c1_RePr
   End Do

   Do i = 1, Nx-1, 1
     c1_dxhmx(i) = 1.0d0 / ( dxhm(i) * dx(i) )
     c1_dxhpx(i) = 1.0d0 / ( dxhp(i) * dx(i) )
   End Do


   !--- y ---
   Do j = 0, Ny+2, 1
     !y(j) = y(j-1) + Ly/Real(Ny+2)!
     !y(j) = Ly * ( 1.0d0 - tanh( alpha_y * (1.0d0 - Dble(j) / Dble(Ny-1)) ) / tanh(alpha_y) )
     y(j) = 0.5d0*Ly*(1.0d0-atanh(tanh(alpha_y)*(1.0d0-2.0d0*(Dble(j)/Dble(Ny-1)))) / alpha_y )
   End Do

   Do j = 1, Ny+2, 1
     dy(j) = y(j) - y(j-1)  
   End Do
     dy(0) = y(0) - ( 2.0d0 * y(0) - y(1) )

   Do j = 0, Ny, 1
     dym(j)    = 0.5d0 * dy(j)
     dyp(j)    = 0.5d0 * dy(j)

     c1_dy (j) = 1.0d0 /  dy(j)
     c1_dy2(j) = 1.0d0 / (dy(j  ) * dy(j))

     c1_dym(j) = 1.0d0 /   dym (j)
     c1_dyp(j) = 1.0d0 /   dyp (j)

   End Do

   Do j = 1, Ny-1, 1
     dyhm(j)     = ( dy(j-1) + dy(j  ) ) * 0.5d0
     dyhp(j)     = ( dy(j  ) + dy(j+1) ) * 0.5d0

     c1_dyhm(j)  = 1.0d0 /   dyhm(j)
     c1_dyhp(j)  = 1.0d0 /   dyhp(j)
   End Do

   Do j = 1, Ny-1, 1
     c1_dyRe(j)   = c1_dy  (j) * c1_Re
     c1_dyhpRe(j) = c1_dyhp(j) * c1_Re
     c1_dyRePr(j) = c1_dy  (j) * c1_RePr
   End Do

   Do j = 1, Ny-1, 1
     c1_dyhmy(j) = 1.0d0 / ( dyhm(j) * dy(j) )
     c1_dyhpy(j) = 1.0d0 / ( dyhp(j) * dy(j) )
   End Do


   !--- z ---
   Do k = 0, Nz+2, 1
     z(k) = z(k-1) + Lz/Real(Nz+2)!
   End Do

   dz = Lz / Dble( Nz-1 )


   !--- Parameters ---
   c1_pi     = 1.0d0 / pi

   c1_Re     = 1.0d0 /  Re
   c1_RePr   = 1.0d0 / (Re * Pr)

   c1_Nx     = 1.0d0 / ( Dble(Nx-1) )                
   c1_Ny     = 1.0d0 / ( Dble(Ny-1) )                
   c1_Nz     = 1.0d0 / ( Dble(Nz-1) )                
   c1_NxNy   = 1.0d0 / ( Dble(Nx-1)*Dble(Ny-1) ) 
   c1_NxNz   = 1.0d0 / ( Dble(Nx-1)*Dble(Nz-1) ) 

   c1_dz     = 1.0d0 /   dz
   c1_dz2    = 1.0d0 /  (dz * dz)
   c1_dzRe   = 1.0d0 /  (dz * Re)
   c1_dz2Re  = 1.0d0 /  (dz * dz * Re)

   c1_dzRePr = 1.0d0 /  (dz *      Re * Pr)
   c1_dz2RePr= 1.0d0 /  (dz * dz * Re * Pr)

   c1_Fr2    = 1.0d0 * c1_Fr * c1_Fr


   Return
End Subroutine


!=======================================================================================================================
!========= [ Subroutine Set_IBM ] ======================================================================================
!=======================================================================================================================
Subroutine Set_IBM

   Use Definition
   Implicit none 
   Integer i, j, k


   !=== Coordinate of U,V,W,Ca ===
   !--- U ---
   X_U(:) = 0.0d0
   Y_U(:) = 0.0d0
   Z_U(:) = 0.0d0
   Do i=0,Nx+2
     x_U(i) = x(i)
   End Do

   Do j=1,Ny+2
     y_U(j) = y(j-1) + c1_2*( y(j)-y(j-1) )
   End Do
   y_U(0) = y_U(1) - ( y_U(2)-y_U(1) )

   Do k=-1,Nz+2
     z_U(k) = z(k) + c1_2*( z(k)-z(k-1) )
   End Do
   z_U(0) = z_U(1) + ( z_U(2)-z_U(1) )


   !--- V ---
   X_V(:) = 0.0d0
   Y_V(:) = 0.0d0
   Z_V(:) = 0.0d0
   Do j=0,Ny+2
     y_V(j) = y(j)
   End Do

   Do i=1,Nx+2
     x_V(i) = x(i-1) + c1_2*( x(i)-x(i-1) )
   End Do
   x_V(0) = x_V(1) - ( x_V(2)-x_V(1) )

   Do k=-1,Nz+2
     z_V(k) = z(k) + c1_2*( z(k)-z(k-1) )
   End Do
   z_V(0) = z_V(1) + ( z_V(2)-z_V(1) )


   !--- W ---
   X_W(:) = 0.0d0
   Y_W(:) = 0.0d0
   Z_W(:) = 0.0d0
   Do k=-2,Nz+2
     z_W(k) = z(k)
   End Do

   Do i=0,Nx+2
     x_W(i) = x_V(i)
   End Do

   Do j=0,Ny+2
     y_W(j) = y_U(j)
   End Do


   !--- Ca ---
   X_Ca(:) = 0.0d0
   Y_Ca(:) = 0.0d0
   Z_Ca(:) = 0.0d0
   Do i=0,Nx+2
     x_Ca(i) = x_V(i)
   End Do

   Do j=0,Ny+2
     y_Ca(j) = y_U(j)
   End Do

   Do k=-2,Nz+2
     z_Ca(k) = z_U(k)
   End Do


   Call Generate_Object

   !=== U-Mark Boundary Grid ===
   Call Mark_BG(y_U ,y_ObUU ,y_ObDU ,j_ObUU ,j_ObDU ,IB_U ,FS_U )
   Call Mark_BG(y_V ,y_ObUV ,y_ObDV ,j_ObUV ,j_ObDV ,IB_V ,FS_V )
   Call Mark_BG(y_W ,y_ObUW ,y_ObDW ,j_ObUW ,j_ObDW ,IB_W ,FS_W )
   Call Mark_BG(y_Ca,y_ObUCa,y_ObDCa,j_ObUCa,j_ObDCa,IB_Ca,FS_Ca)

   !=== Set Imaging Point ===
   Call Set_IP(x_ObOU,y_ObOUU,y_ObODU,x_IPU,y_IPUU,y_IPDU,x_U,y_U,y_ObUU,y_ObDU,j_ObUU,j_ObDU )


   !--- Grid ---
   Call Write_Grid(1,x_ObU, y_ObUU ,y_ObDU ,x_U ,y_U )
   Call Write_Grid(2,x_ObV ,y_ObUV ,y_ObDV ,x_V ,y_V )
   Call Write_Grid(3,x_ObW ,y_ObUW ,y_ObDW ,x_W ,y_W )
   Call Write_Grid(4,x_ObCa,y_ObUCa,y_ObDCa,x_Ca,y_Ca)

   !--- Mark_Boundary Grid ---
   Call Write_MarkBG(1,j_ObUU ,j_ObDU ,x_U ,y_U )
   Call Write_MarkBG(2,j_ObUV ,j_ObDV ,x_V ,y_V )
   Call Write_MarkBG(3,j_ObUW ,j_ObDW ,x_W ,y_W )
   Call Write_MarkBG(4,j_ObUCa,j_ObDCa,x_Ca,y_Ca)

   !--- Imaging Points ---


   Return
End Subroutine


!=======================================================================================================================
!========= [ Subroutine Generate_Object ] ==============================================================================
!=======================================================================================================================
Subroutine Generate_Object

   Use Definition
   Implicit none 
   Integer i, j, k

   Real*4  ObShape_U, ObShape_D

   !--- Object Range ---
   Do i=0,Nx-1
     If (x(i).GE.(x_Ob0-Lx_Ob)) Then
       i_OL = i
       Exit
     End If
   End Do

   Do i=i_OL,Nx-1
     If (x(i).GE.(x_Ob0+Lx_Ob)) Then
       i_OR = i
       Exit
     End If
   End Do

   Do j=0,Ny-1
     If (y(j).GE.(y_Ob0-Ly_Ob)) Then
       j_OD = j-1
       Exit
     End If
   End Do

   Do j=j_OU,Ny-1
     If (y(j).GE.(y_Ob0+Ly_Ob)) Then
       j_OU = j
       Exit
     End If
   End Do


   !--- Generate Object ---
   Do i=i_OL,i_OR
     x_ObU(i)  = x_U(i)
     x_ObV(i)  = x_V(i)
     x_ObW(i)  = x_W(i)
     x_ObCa(i) = x_Ca(i)

     y_ObUU(i)  = ObShape_U(x_ObU ,i)
     y_ObUV(i)  = ObShape_U(x_ObV ,i)
     y_ObUW(i)  = ObShape_U(x_ObW ,i)
     y_ObUCa(i) = ObShape_U(x_ObCa,i)

     y_ObDU(i)  = ObShape_D(x_ObU ,i)
     y_ObDV(i)  = ObShape_D(x_ObV ,i)
     y_ObDW(i)  = ObShape_D(x_ObW ,i)
     y_ObDCa(i) = ObShape_D(x_ObCa,i)
   End Do



   Return
End Subroutine


!==========================================================================================================
!========= [ Subroutine Mark_BG ] =========================================================================
!==========================================================================================================
Subroutine Mark_BG(y1,yu,yd,ju,jd,IB,FS)

 Use Definition
 Implicit none

 Integer i,j,jj,k,ju(-2:Nx+2),jd(-2:Nx+2),IB(-2:Nx+2,-2:Ny+2),FS(-2:Nx+2,-2:Ny+2)
 Integer j1,j2
 Real*4  y1(-2:Ny+2),yu(-2:Nx+2),yd(-2:Nx+2)


 !=== 固体内の格子点を探索 ===
 !--- 境界に最も近傍に位置する固体内の格子点に印付け ---
 ju(:) = 0
 jd(:) = 0
 Do i=i_OL,i_OR
 Do j=0,Ny-1
   If (y1(j).LT.yu(i)) ju(i) = j	!iに対する境界近傍の点のj
 End DO ; End Do

 Do i=i_OL,i_OR
   Do j=0,Ny-1
     If (y1(j).GT.yd(i)) Then
       jd(i) = j	!iに対する境界近傍の点のj
       Exit
     End If
   End DO
 End Do


 !--- 流体領域と固体領域の印付け ---
 IB(:,:) = 1 !初期値:全流体領域
 FS(:,:) = 1 !初期値:全流体領域

 Do i=i_OL,i_OR
   j1 = jd(i)
   j2 = ju(i)

   Do j=j1,j2
     IB(i,j) = 0 !固体領域
   End Do

 End Do



 Return
End Subroutine Mark_BG


!==========================================================================================================
!========= [ Subroutine Set_IP ] ==========================================================================
!==========================================================================================================
Subroutine Set_IP(xobo,yobUo,yobDo,xp,ypU,ypD,x1,y1,yobU,yobD,jobU,jobD)
 
   Use Definition
   Implicit none

   Real*4 xobo(-2:Nx+2),yobUo(-2:Nx+2),yobDo(-2:Nx+2),xp(-2:Nx+2),ypU(-2:Nx+2),ypD(-2:Nx+2)
   Real*4 x1(-2:Nx+2),y1(-2:Ny+2),yobU(-2:Nx+2),yobD(-2:Nx+2)
   Real*4 xmo(-2:Nx+2),ymo(-2:Nx+2)
   Real*4 ObShape_U,ObShape_D
   Real*4 ObShape_UIPU,ObShape_UIPD
   Integer jobU(-2:Nx+2),jobD(-2:Nx+2),i,j,k

   Real*4  s,f,h,fd,dms0
   Real*4 ,Parameter :: ds = 1.0d-3


   !--- 固体内格子と境界が直交するときの境界の座標 ---
   xobo (:) = 0.0d0
   yobUo(:) = 0.0d0
   yobDo(:) = 0.0d0
   j = 0

   Do i=i_OL,i_OR

     !上側--------------------------------------------
     j = jobU(i)
     s = x(i)
     Do k=1,100000

 	f = ObShape_UIPU(s   ,y1,j,x1,i)
 	h = ObShape_UIPU(s+ds,y1,j,x1,i)

 	fd=(h-f)/ds

 	dms0 = -f / fd
 	   s = s + dms0

	if(abs(dms0/s).LT.1.d-6) Exit
     End Do

     xobo (i) = s
     yobDo(i) = ObShape_U(xobo,i)!0.5d0*( 1.0d0+cos(xmo(i)-x(i_ML)+PI) )
     !--------------------------------------------


     !下側--------------------------------------------
     j = jobD(i)
     s = x(i)
     Do k=1,100000

 	f = ObShape_UIPD(s   ,y1,j,x1,i)
 	h = ObShape_UIPD(s+ds,y1,j,x1,i)

 	fd=(h-f)/ds

 	dms0 = -f / fd
 	   s = s + dms0

	if(abs(dms0/s).LT.1.d-6) Exit
     End Do

     xobo (i) = s
     yobUo(i) = ObShape_D(xobo,i)!0.5d0*( 1.0d0+cos(xmo(i)-x(i_ML)+PI) )
     !--------------------------------------------

   End Do





   Return
End Subroutine Set_IP


!==========================================================================================================
!========= [ Subroutine Write_Grid ] ======================================================================
!==========================================================================================================
Subroutine Write_Grid(n,xob,yobU,yobD,x1,y1)
 
   Use Definition
   Implicit none

   Integer i,j,k,n
   Real*4  xob(-2:Nx+2),yobU(-2:Nx+2),yobD(-2:Nx+2),x1(-2:Nx+2),y1(-2:Ny+2)


   Character filename*128

   Write(filename, '("GR_Gene/12_Object", i2.2, ".txt")') n
   Open (12,file=filename)

   Write(filename, '("GR_Gene/13_Velocity_A", i2.2, ".txt")') n
   Open (13,file=filename)

   Write(filename, '("GR_Gene/14_Velocity_A", i2.2, ".txt")') n
   Open (14,file=filename)


   !--- Mountain ---
   Do i=i_OL,i_OR
     Write(12,*) xob(i),yobU(i)
   End Do
   !---

   Write(12,*)

   Do i=i_OL,i_OR
     Write(12,*) xob(i),yobD(i)
   End Do


   !--- U1 ---
   Do j=0,Ny-1
   Do i=0,Nx-1
     Write(13,*) x1(i),y1(j)!,z1(k)
     If (i.EQ.Nx-1) Write(13,*)
   End Do ; End Do
   !---

   !--- U2 ---
   Do i=0,Nx-1
   Do j=0,Ny-1
     Write(14,*) x1(i),y1(j)
     If (j.EQ.Ny-1) Write(14,*)
   End Do ; End Do
   !---




   Return
End Subroutine


!==========================================================================================================
!========= [ Subroutine Write_MarkBG ] ====================================================================
!==========================================================================================================
Subroutine Write_MarkBG(n,j_obU,j_obD,x1,y1)

   Use Definition
   Implicit none

   Integer n, j_obU(-2:Nx+2), j_obD(-2:Nx+2)
   Integer i,j,k,j1,j2
   Real*4  x1(-2:Nx+2),y1(-2:Ny+2)

   Character filename*128

   Write(filename, '("GR_Gene/15_MarkBG", i2.2, ".txt")') n
   Open (15,file=filename)

   !--- Mountain & Boundary Grid ---
   Do i=i_OL,i_OR
     j1 = j_ObU(i)
     j2 = j_ObD(i)
     Write(15,*) x1(i),y1(j1),y1(j2)
   End Do


   Return
End Subroutine


!==========================================================================================================
!========= [ Function ObShape_U ] =========================================================================
!==========================================================================================================
Function ObShape_U(x1,i)

   Use Definition
   Implicit none
   Integer i,j,k,n

   Real*4 x1(-2:Nx+2), ObShape_U, xx

   xx = Ly_Ob**2.0d0 - (x1(i)-x_Ob0)**2.0d0
   ObShape_U = ( dmax1( 0.00001d0,xx ) )**0.5d0 + y_Ob0

   If (i.EQ.i_OL.OR.i.EQ.i_OR) ObShape_U = y_Ob0 + eps

   Return
End Function


!==========================================================================================================
!========= [ Function ObShape_U ] =========================================================================
!==========================================================================================================
Function ObShape_D(x1,i)

   Use Definition
   Implicit none
   Integer i,j,k,n

   Real*4 x1(-2:Nx+2), ObShape_D, xx


   xx = Ly_Ob**2.0d0 - (x1(i)-x_Ob0)**2.0d0
   ObShape_D = - ( dmax1( 0.00001d0,xx ) )**0.5d0 + y_Ob0

   If (i.EQ.i_OL.OR.i.EQ.i_OR) ObShape_D = y_Ob0 - eps

   Return
End Function


!==========================================================================================================
!========= [ Function ObShape_U ] =========================================================================
!==========================================================================================================
Function ObShape_UIPU(s,y1,j,x1,i)

   Use Definition
   Implicit none
   Integer i,j,k,n

   Real*4 x1(-2:Nx+2), y1(-2:Ny+2), ObShape_UIPU, xx, s


   xx = Ly_Ob**2.0d0 - (s-x_Ob0)**2.0d0
   xx = ( dmax1( 0.00001d0,xx ) )**0.5d0
   ObShape_UIPU = ( xx ) * ( y1(j) - ( xx+y_Ob0 ) ) + ( x1(i)-s )


   Return
End Function


!==========================================================================================================
!========= [ Function ObShape_U ] =========================================================================
!==========================================================================================================
Function ObShape_UIPD(s,y1,j,x1,i)

   Use Definition
   Implicit none
   Integer i,j,k,n

   Real*4 x1(-2:Nx+2), y1(-2:Ny+2), ObShape_UIPD, xx, s


   ObShape_UIPD = 1.0d0


   Return
End Function
