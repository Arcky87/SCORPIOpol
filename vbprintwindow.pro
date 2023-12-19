PRO VBPrintWindow, DrawId
     ;             .
       ;           .
         ;         .
  ;Get the window index of the drawable to be printed
  WIDGET_CONTROL, DrawId, Get_Value=Index
    ;             .
       ;           .
        ;          .
  ;Create a Printer object and draw the graphic to it
  oPrinter = OBJ_NEW ('IDLgrPrinter')

  ;Display a print dialog box
  Result = DIALOG_PRINTERSETUP(oPrinter)
       ;           .
           ;       .
             ;     .
  oPrinter->Draw, oView
            ;      .
              ;    .
               ;   .
END ;VBPrintWindow

