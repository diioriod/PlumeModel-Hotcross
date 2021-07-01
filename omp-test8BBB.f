
      program test8bbb

! pgf95 -mp -r8 -O3 -Mextend omp-test8BBB.f

      INTEGER n_threads, tid, OMP_get_num_threads, OMP_get_thread_num

      integer nx,ny,nz

      parameter (nx = 256, ny= 256, nz = 140)

      real, allocatable:: a(:,:,:), b(:,:,:), c(:,:,:) 
      real time1, time2 
      logical answer

      allocate ( a(0:nx+1,0:ny+1,0:nz+1),b(0:nx+1,0:ny+1,0:nz+1),c(0:nx+1,0:ny+1,0:nz+1))

!LOAD DUMMY DATA

         a = 5; b= 6; c= 7
      do k = 1, nz
        do j = 1, ny
         do i = 1, nx
           a(i,j,k) = i/real(j)*k; b(i,j,k) = 6*k*sqrt(real(i))/real(j) ; c(i,j,k)=7*i*(k*j)**1.1d0
      enddo;enddo;enddo

!FINISH LOAD OF DUMMY DATA

      write(6,*) 'number of threads    execution time(s) '

      do nthreads = 1, 16

         time1 = omp_get_wtime() 
         call OMP_set_num_threads(nthreads)

         mem = 1
         call sedona (a,b,c,nx,ny,nz,mem)

         mem = 2
         call sedona (b,a,c,nx,ny,nz,mem)

         time2 = omp_get_wtime()

!$OMP PARALLEL

         if ( OMP_get_thread_num() .eq. 0) then
            write(6,*) omp_get_num_threads(), time2-time1
         endif

!$OMP END PARALLEL

      enddo

      end program test8bbb

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine sedona(A,B,C,nx,ny,nz,mem)
      integer nx,ny,nz,mem

      real, intent(in ):: a(0:nx+1,0:ny+1,0:nz+1),b(0:nx+1,0:ny+1,0:nz+1)
      real, intent(out):: c(0:nx+1,0:ny+1,0:nz+1)      
      integer i,j,k, tid,  OMP_get_thread_num
      real time 

!$OMP PARALLEL  private(tid) Shared(a,b,c) 

!$OMP DO 

       do l = 1, 100
       do k = 1, nz
        do j = 1, ny
         do i = 1, nx 
               c(i,j,k) =  a(i-1,j,k) * A(i,j,k) + A(i+1,j,k)* A(i,j-1,k)
     &                  + A(i,j+1,k)+ sin(A(i,j,k+1))+ A(i,j,k-1) + l
       enddo;enddo;enddo;enddo
       
!$OMP END DO 

!$OMP END PARALLEL 


      end  subroutine sedona