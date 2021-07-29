program Ising2D
  
    implicit none
    integer :: nsite,i,k,im,ip,km,kp
    parameter(nsite = 10)
    real :: ispin(nsite,nsite),xx,ifield
    real :: AJ=1.0,B=0.0
    real :: iU0=0.0,U0=0.0
  
   ! 初期スピン配置の生成
   
    do i=1,nsite
       do k=1,nsite
          ispin(i,k)=1
       enddo
    enddo

    do i=1,nsite
        do k=1,nsite
          call random_number(xx)
          if (xx > 0.5d0 ) then
            ispin(i,k) = -1
          endif
 
          kp=k+1
          if(kp == nsite+1) kp=1
          km=k-1
          if(km == 0) km=nsite
          ip=i+1
          if(ip == nsite+1) ip=1
          im=i-1
          if(im == 0) im=nsite
          
    
          ifield=ispin(ip,k)+ispin(im,k)+ispin(i,kp)+ispin(i,km) !ΣSiの計算
       
          iU0 = (-AJ * ifield - B) * ispin(i,k)     ! 各サイトごとのU0を計算
          U0 = U0 + iU0
          write (*,*) iU0,U0
       end do
    end do
end program
