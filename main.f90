
      program main

      implicit none

      real*8 :: x0(3,4),x1(3,4)
      integer :: ifn(3,4) 
      integer :: n,j
      real*8 :: dV,r

      ifn(:,1) = (/1, 4, 3/)
      ifn(:,2) = (/1, 3, 2/)
      ifn(:,3) = (/2, 3, 4/)
      ifn(:,4) = (/1, 2, 4/)

      x0(:,1) = (/1.0d0, 0.0d0, 0.0d0/)
      x0(:,2) = (/0.0d0, 1.0d0, 0.0d0/)
      x0(:,3) = (/0.0d0, 0.0d0, 0.0d0/)
      x0(:,4) = (/0.0d0, 0.0d0, 1.0d0/)

      x1 = x0
      do n=1,4
         call perturb_node(x1(:,n))
      enddo

      dV = tet_vol(x1(:,1),x1(:,2),x1(:,3),x1(:,4))
      dV = dV - tet_vol(x0(:,1),x0(:,2),x0(:,3),x0(:,4))

      PRINT*,'volume change = ',dV

      dV = 0.0d0
      do j=1,4
         dV = dV + tri_dV(x0(:,ifn(:,j)),x1(:,ifn(:,j)))
      enddo
      
      PRINT*,'volume change = ',dV
 
      contains

      subroutine perturb_node(x)
      implicit none
      real(8), intent(INOUT) :: x(3)
      real(8) :: r
      integer :: n

      do n=1,3
         call random_number(r)
         r = r - 0.5 ! between -0.5 and 0.5
         r = r*0.5d0 ! between -0.25 and 0.25
         x(n) = x(n) + r
      enddo
        
      end subroutine 

      function tri_dV(x0,x1) result(dV)
      implicit none
      real(8), intent(IN) :: x0(3,3),x1(3,3)
      real*8 :: w0(3),w1(3),w2(3),w3(3)
      real(8) :: Savg(3),xavg(3,3),dV
      integer :: i

      w1(:) = x1(:,1) - x0(:,1)
      w2(:) = x1(:,2) - x0(:,2)
      w3(:) = x1(:,3) - x0(:,3)

      w0 = (w1 + w2 + w3)/3.0d0
      do i=1,3
          xavg(:,i) = 0.5d0*(x0(:,i) + x1(:,i))
      enddo

      Savg(:) = tri_face( xavg(:,1), xavg(:,2), xavg(:,3) )

      dV = dot_product(w0,Savg) + dot_product(w1,cross_product(w2,w3))/24.0d0

      end function

      function cross_product(v1,v2) result(v3)
      implicit none
      real(8), intent(IN) :: v1(3),v2(3)
      real(8) :: v3(3)
      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      end function

      function tet_vol(x1,x2,x3,x4) result(vol)
      implicit none
      real(8), intent(IN) :: x1(3),x2(3),x3(3),x4(3)
      real(8) :: vol
      real(8), dimension(3) :: v41,v42,v43,c12
      v41(:) = x1(:) - x4(:)    
      v42(:) = x2(:) - x4(:)  
      v43(:) = x3(:) - x4(:)
      c12 = cross_product(v42,v41) 
      vol = dot_product(c12,v43)/6.0d0 
      end function 

      function tri_face(x1,x2,x3) result(face)
      implicit none
      real(8), intent(IN) :: x1(3),x2(3),x3(3)
      real(8) :: c1(3),v1(3),v2(3),face(3)
      v1 = x2(:) - x1(:)
      v2 = x3(:) - x2(:)
      face = 0.5d0*cross_product(v1,v2)
      return
      end function

      end program
