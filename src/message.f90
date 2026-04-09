subroutine error_message (loc, text)
    character(*) :: loc, text
    write (*,*) "ERROR (", loc, "): ", text
    stop
end subroutine error_message