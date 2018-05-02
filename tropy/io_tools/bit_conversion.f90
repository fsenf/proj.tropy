      !--------------------------------------------------------------------!
 
      subroutine bit_conversion(dat8bit, N8bit, dat10bit)
      implicit none
 
      !f2py integer(KIND=4),INTENT(IN) :: N8bit
      !f2py integer(KIND=4), DIMENSION(N8bit),INTENT(IN) :: dat8bit
      !f2py integer(KIND=4), DIMENSION(8 * N8bit / 10),INTENT(OUT) :: dat10bit

  
      integer(KIND=4),INTENT(IN) :: N8bit
      integer(KIND=4), DIMENSION(N8bit),INTENT(IN) :: dat8bit
      integer(KIND=4), DIMENSION(8 * N8bit / 10),INTENT(OUT) :: dat10bit

      integer(KIND=4) :: Nbits, N10bit
      integer(KIND=4) :: n, n8, n10, i, res
 
      integer(KIND=4), DIMENSION(8 * N8bit) :: bits

      integer(KIND=4), DIMENSION(8 * N8bit / 10) :: base


      Nbits = N8bit * 8
      N10bit = (N8bit * 8)/ 10

      do i = 1, 10
        base(i) = 2**(i - 1)
      enddo

      n = 0
      do n8 = 1, N8bit
        res = dat8bit(n8)

        do i = 8, 1, -1
          n = n + 1
          
          ! get the bit belonging to base(i) 
          bits(n) = res / base(i)

          ! calculate the residual
          res = res - bits(n) * base(i)
        enddo   
      enddo

 
      dat10bit = 0
      n = 0
      do n10 = 1, N10bit
        do i = 10, 1, -1 
          n = n + 1
          ! sum all up again with base10
          dat10bit(n10) = dat10bit(n10) + bits(n) * base(i)
        enddo
      enddo

      return
      end 
