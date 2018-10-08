        character(len=*) function removeblanks(string)

! Remove embedded blanks from string.

        implicit real (a-h,o-z)
        implicit integer (i,n)
        character(len=*) string

        removeblanks = ''
        l = len(string)
        m = 0

        do i = 1,l
            if (index(string(i:i),' ') == 0) then
                m = m+1
                removeblanks(m:m) = string(i:i)
            endif
        enddo

        return

        end

