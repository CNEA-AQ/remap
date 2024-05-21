module remap_discrete_mod

 !use netcdf 
 use PROJ
 implicit none

 !Remap method ID:
 !interpolation (to finer grid)
 integer :: bilinear          =101 !cont. fields (real/double)
 integer :: bicubic           =102 !cont. fields
 integer :: dist_weighted     =102 !cont. fields
 integer :: nearest_neighboor =102 !categorical  (int)
 !remap (to coarser grid)
 integer :: conserv_1order    =201 !cont. fields (real/double)
 integer :: conserv_2order    =202 !cont. fields
 integer :: average           =203 !cont. fields
 integer :: mode              =204 !categorical  (int)

 interface remap_discrete_field
   module procedure remap_discrete_1d_field, remap_discrete_2d_field
 end interface remap_discrete_field

contains

subroutine remap_discrete_2d_field(arr1,arr2,g1,g2,method)
   implicit none
   integer                , intent(in)         :: arr1(:,:) !src arrays
   integer                , intent(inout)      :: arr2(:,:) !dst arrays (OUT)
   type(regular_grid_type), intent(in)         :: g1, g2    !src and dst grid
   character(*)           , intent(in)         :: method    !bilinear, bicubic, conservative 
   real (SCRIP_r8), allocatable, dimension(:)  :: grid1_array,grid2_array

   grid1_array=reshape(arr1,[g1%nx*g1%ny])
   grid2_array=reshape(arr2,[g2%nx*g2%ny])

   call remap_discrete_1d_field(grid1_array,grid2_array,g1,g2,method)

   arr2=reshape(grid2_array,[g2%nx,g2%ny])
end subroutine

subroutine remap_discrete_1d_field(grid1_array,grid2_array,g1,g2,method)
   implicit none
   integer                , intent(in)         :: arr1(:,:) !src arrays
   integer                , intent(inout)      :: arr2(:,:) !dst arrays (OUT)
   type(regular_grid_type), intent(in)         :: g1, g2    !src and dst grid
   character(*)           , intent(in)         :: method    !bilinear, bicubic, conservative 

end subroutine

function interpolate(p,g,inp_file,varname,method)       result(img2)
 implicit none
 type(grid_type), intent(in)  :: g  !desired grid
 type(proj_type), intent(in)  :: p  !proj of desired grid
 character(*), intent(in)     :: inp_file,varname,method
 integer                      :: methodId
 real,allocatable :: img2(:,:)      !output array

 integer :: i,j,k

 type(grid_type)  :: GG,GC          !global grid (input grid) &  global grid (CROPPED)
 real, allocatable :: img1(:,:)      !cropped img to be interpolated

 integer :: is,ie,js,je     !indices that defines subarray
 real    :: px,py,x,y       !dummy variables for temporal coordinates
 real    :: w11,w12,w21,w22 !weights for bilinear interpolation
 real    :: p11,p12,p21,p22 !params. for bilinear interpolation
 real    :: x1,x2,y1,y2     !params. for interpolation
 integer :: i1,i2,j1,j2     !dummy indexes for interpolation

 integer :: scale_x,scale_y

 print*,"  Interpolando: "//trim(inp_file)//":"//trim(varname)//"..."
 ! Asumo que estoy trabajando con grillas regulares (dx/dy =cte.).    
 ! Asumo que lat y lon estan ordenados de forma creciente.

 call get_cropped_img(p,g,inp_file,varname,img1,GC)
 
 !Veo si la grilla destino es mas densa o no que la original.
 call xy2ll(p,g%xmin,g%ymin,x1,y1)  !
 call xy2ll(p,g%xmax,g%ymax,x2,y2)  !

 scale_x=CEILING((x2-x1)/(g%nx)/GC%dx)
 scale_y=CEILING((y2-y1)/(g%ny)/GC%dy)

 if ( method == "mode"    ) then
     methodId=1
 else if ( method == "near-neighbor"     ) then
     methodId=2

 endif

 allocate(img2(g%nx,g%ny))  !array a interpolar:

 !REGRIDDING:
 if( methodId == 1) then

    do i=1,g%nx
       do j=1,g%ny
           idx= (j-1)*g%nx + (i-1)  !index mapping: 2D -> 1D

           px=g%xmin+g%dx*i  !projected coordinate-x
           py=g%ymin+g%dy*j  !projected coordinate-y
                                                                                       
           call xy2ll(p,px,py,x,y)  !I want coordinates in same proj than global file.

           i1=MAX(     1, FLOOR( (x-GC%lonmin) / GC%dx - scale_x*0.5 ) )
           j1=MAX(     1, FLOOR( (y-GC%latmin) / GC%dy - scale_y*0.5 ) )
           i2=MIN( GC%nx, i1+scale_x                                   )
           j2=MIN( GC%ny, j1+scale_y                                   )

           if ( i1 >= 1 .and. i2 <= GC%nx .and. j1 >=1 .and. j2 <= GC%ny ) then
               img2(idx)=MODE(INT(img1(i1:i2,j1:j2)))                        !mode
           else
               img2(i,j)=0
           endif
       enddo
    enddo
 endif
 !NEAREST-NEIGBOR
 if ( methodId == 2 ) then
 
   do i=1,g%nx
      do j=1,g%ny
         idx= (j-1)*g%nx + (i-1)  !index mapping: 2D -> 1D

         px=g%xmin+g%dx*i  !projected coordinate-x                                    
         py=g%ymin+g%dy*j  !projected coordinate-y
                                                                                      
         call xy2ll(p,px,py,x,y)  !I want coordinates in same proj than global file.

        ii=floor((x - cg%lonmin) / gc%dx)
        jj=floor((y - cg%latmin) / gc%dy)
        img2(i,j)=img1(ii,jj)
         
        iidx= (jj-1)*gc%nx + (ii-1)
        img2(idx)=img1(iidx)

      enddo
   enddo

 endif

end function


end module

