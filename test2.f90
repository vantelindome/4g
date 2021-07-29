program rand
implicit none
integer :: nsite,i,k,im,ip,km,kp
parameter(nsite=10)
real(4) :: xx
integer :: try
real :: ispin(nsite,nsite)
real :: AJ=1.0,B=0.0
integer :: ltemp,T=2
real :: H,bratio,deltaU, U, totalU, total2U
integer :: x
real :: C,kB=1.0


   do ltemp=1,T
        deltaU = 0
        U = 0
        totalU = 0
        total2U = 0

        do try=1,200
            call random_number(xx)
            if (xx >= 0.0 .and. xx < 0.1) x=1
            if (xx >= 0.1 .and. xx < 0.2) x=2
            if (xx >= 0.2 .and. xx < 0.3) x=3
            if (xx >= 0.3 .and. xx < 0.4) x=4
            if (xx >= 0.4 .and. xx < 0.5) x=5
            if (xx >= 0.5 .and. xx < 0.6) x=6
            if (xx >= 0.6 .and. xx < 0.7) x=7
            if (xx >= 0.7 .and. xx < 0.8) x=8
            if (xx >= 0.8 .and. xx < 0.9) x=9
            if (xx >= 0.9 .and. xx < 1.0) x=10
                
            i=x

            call random_number(xx)
            if (xx >= 0.0 .and. xx < 0.1) x=1
            if (xx >= 0.1 .and. xx < 0.2) x=2
            if (xx >= 0.2 .and. xx < 0.3) x=3
            if (xx >= 0.3 .and. xx < 0.4) x=4
            if (xx >= 0.4 .and. xx < 0.5) x=5
            if (xx >= 0.5 .and. xx < 0.6) x=6
            if (xx >= 0.6 .and. xx < 0.7) x=7
            if (xx >= 0.7 .and. xx < 0.8) x=8
            if (xx >= 0.8 .and. xx < 0.9) x=9
            if (xx >= 0.9 .and. xx < 1.0) x=10
                
            k=x

            !周期的境界条件
            kp=k+1
            if(kp == nsite+1) kp=1
            km=k-1
            if(km == 0) km=nsite
            ip=i+1
            if(ip == nsite+1) ip=1
            im=i-1
            if(im == 0) im=nsite
                    
            !メトロポリス判定     
            if(deltaU <= 0.0d0) then
                ispin(i,k) = (-1)*ispin(i,k)
            else
                call random_number(xx)
                bratio = exp(-deltaU/ltemp)
                if(xx < bratio) ispin(i,k) = (-1)*ispin(i,k)
            endif   

            if(try>100) then
                H = ispin(ip,k) + ispin(im,k) + ispin(i,kp) + ispin(i,km) !ΔUの計算
            
                deltaU = (-AJ * H - B)*ispin(i,k)            ! Eの計算
                U = U + deltaU
                totalU = totalU + U
                total2U = total2U + U * U 
                write(*,*) deltaU,U,totalU,total2U
            endif

          
        enddo   

        C = (total2U - (totalU * totalU)) / (kB * ltemp * ltemp)
        write(*,*) ltemp, totalU, C
    enddo

end program
