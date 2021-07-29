program Ising2D
    implicit none
    integer :: nsite,i,k,im,ip,km,kp
    parameter(nsite = 20)
    integer :: TT
    parameter(TT=300000)
    real :: ispin(nsite,nsite),xx,ifield
    real :: AJ=1.0,B=0.0
    real :: iU0=0.0
    integer :: ltemp,T=400
    integer :: try
    real :: Ei,Ek
    real :: H,deltaU, U=0.0, totalU=0.0, total2U=0.0
    integer :: x
    real :: C=0,kB=1.0
    real :: bratio
            
       ! 初期スピン配置の生成
             
     do i=1,nsite
        do k=1,nsite
           ispin(i,k)=1            !すべてのスピンを1とする
        enddo
     enddo
        
     do i=1,nsite
         do k=1,nsite
           call random_number(xx)
           if (xx > 0.5d0 ) then
              ispin(i,k) = -1      !確率1/2でスピンを―1にする
           endif
    
           !周期的境界条件
           kp=k+1
           if(kp == nsite+1) kp=1
           km=k-1
           if(km == 0) km=nsite
           ip=i+1
           if(ip == nsite+1) ip=1
           im=i-1
           if(im == 0) im=nsite
             
           ifield = ispin(ip,k) + ispin(im,k) + ispin(i,kp) + ispin(i,km) !ΣSiの計算
                 
           iU0 = (-AJ * ifield -B) * ispin(i,k)     ! 各サイトごとのUを計算
           U = U + iU0                       ! 初期エネルギーUができる
        end do
     end do
           
     !各温度ごとにMCシミュレーションを回す
     do ltemp=1,T
        do try=1,TT
           call random_number(xx)
           if (xx >= 0.0 .and. xx < 0.05) x=1
           if (xx >= 0.05 .and. xx < 0.1) x=2
           if (xx >= 0.1 .and. xx < 0.15) x=3
           if (xx >= 0.15 .and. xx < 0.2) x=4
           if (xx >= 0.2 .and. xx < 0.25) x=5
           if (xx >= 0.25 .and. xx < 0.3) x=6
           if (xx >= 0.3 .and. xx < 0.35) x=7
           if (xx >= 0.35 .and. xx < 0.4) x=8
           if (xx >= 0.4 .and. xx < 0.45) x=9
           if (xx >= 0.45 .and. xx < 0.5) x=10
           if (xx >= 0.5 .and. xx < 0.55) x=11
           if (xx >= 0.55 .and. xx < 0.6) x=12
           if (xx >= 0.6 .and. xx < 0.65) x=13
           if (xx >= 0.65 .and. xx < 0.7) x=14
           if (xx >= 0.7 .and. xx < 0.75) x=15
           if (xx >= 0.75 .and. xx < 0.8) x=16
           if (xx >= 0.8 .and. xx < 0.85) x=17
           if (xx >= 0.85 .and. xx < 0.9) x=18
           if (xx >= 0.9 .and. xx < 0.95) x=19
           if (xx >= 0.95 .and. xx < 1.0) x=20
                    
           i=x
        
           call random_number(xx)
           if (xx >= 0.0 .and. xx < 0.05) x=1
           if (xx >= 0.05 .and. xx < 0.1) x=2
           if (xx >= 0.1 .and. xx < 0.15) x=3
           if (xx >= 0.15 .and. xx < 0.2) x=4
           if (xx >= 0.2 .and. xx < 0.25) x=5
           if (xx >= 0.25 .and. xx < 0.3) x=6
           if (xx >= 0.3 .and. xx < 0.35) x=7
           if (xx >= 0.35 .and. xx < 0.4) x=8
           if (xx >= 0.4 .and. xx < 0.45) x=9
           if (xx >= 0.45 .and. xx < 0.5) x=10
           if (xx >= 0.5 .and. xx < 0.55) x=11
           if (xx >= 0.55 .and. xx < 0.6) x=12
           if (xx >= 0.6 .and. xx < 0.65) x=13
           if (xx >= 0.65 .and. xx < 0.7) x=14
           if (xx >= 0.7 .and. xx < 0.75) x=15
           if (xx >= 0.75 .and. xx < 0.8) x=16
           if (xx >= 0.8 .and. xx < 0.85) x=17
           if (xx >= 0.85 .and. xx < 0.9) x=18
           if (xx >= 0.9 .and. xx < 0.95) x=19
           if (xx >= 0.95 .and. xx < 1.0) x=20
                    
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
           H = ispin(ip,k) + ispin(im,k) + ispin(i,kp) + ispin(i,km)   ! 最近接スピンの計算
           Ei = (-AJ * H - B) * ispin(i,k)
    
           ispin(i,k) = (-1) * ispin(i,k)           ! 選択した(i,k)についてスピンの正負を変更
           Ek = (-AJ * H - B) * ispin(i,k)          ! 変更した部分のエネルギーを計算
           deltaU =  Ek-Ei                          ! ΔUの計算
    
        !メトロポリス判定     
           if(deltaU <= 0.0d0) then
              ispin(i,k) = 1*ispin(i,k)
           else
              bratio=exp(-deltaU/ltemp)
              call random_number(xx)
              if(xx < bratio) then
                    ispin(i,k) = 1*ispin(i,k)
               else 
                  ispin(i,k) = (-1)*ispin(i,k)
                  deltaU = 0.0
              endif
           endif   
    
           if(try>100000) then
               U = U + deltaU
               totalU = totalU + U
               total2U = total2U + U * U 
             
           endif             
        enddo
        totalU = totalU/200000
        C = (total2U - (totalU * totalU)) / (kB * ltemp * ltemp * 200000 *200000)
        write(*,*) ltemp, totalU, C
 
        deltaU = 0.0
        U = 0.0
        totalU = 0.0
        total2U = 0.0
 
        !終状態を次の温度についての初期状態とする
        do i=1,nsite
            do k=1,nsite
       
              !周期的境界条件
              kp=k+1
              if(kp == nsite+1) kp=1
              km=k-1
              if(km == 0) km=nsite
              ip=i+1
              if(ip == nsite+1) ip=1
              im=i-1
              if(im == 0) im=nsite
                
              ifield = ispin(ip,k) + ispin(im,k) + ispin(i,kp) + ispin(i,km) !ΣSiの計算
                    
              iU0 = (-AJ * ifield -B) * ispin(i,k)     ! 各サイトごとのUを計算
              U = U + iU0                       ! 初期エネルギーUができる
           enddo
     enddo
    enddo
 
 end program Ising2D
 