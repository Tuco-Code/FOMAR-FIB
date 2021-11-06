#include "colors"
#include "textures"



/*{{{ camera, llum, fons */

light_source { <-25000,100,2000>, 1.5}

light_source { <0,5000,2000>, 1.5}

camera {
  location  <5, 10 , 45>
  right -4/3*x
  angle 20
  look_at   < 6.0,1,0>
}
background{rgb .02}

 plane { <0, 1, 0>, -0.5     // normal, displacement
    pigment {
      checker color White, color Gray
    }
  }

/*{{{ Pla inclinat inexistent */
/*
#declare Spla = sin(45*pi/180);
#declare Cpla = cos(45*pi/180);

 plane { <-Spla, Cpla, 0>, -4.0     // normal, displacement
    pigment {
      checker color Red  , color Yellow
    }
  }
*/
/*}}}*/

/*}}}*/

/*{{{ Definicio de sistema de referencia*/
#declare SRef=

   union{
         cylinder{<0, 0, 0>, <1.50, 0, 0>, .025 texture { pigment { color Red   } } }
   //    cylinder{<0, 0, 0>, <0, 1.40, 0>, .025 texture { pigment { color Green } } }
         cylinder{<0, 0, 0>, <0, 0, 1.50>, .025 texture { pigment { color Blue  } } }
        }
        
/*}}}*/


/*{{{ Lectura Taula Denavit-Hartenberg: */

//      alfa    a    |     d    fi      // 

#fopen Dades "taula-DH" read

#declare Taula = array[4][10];
# declare NumElements = 0;

#while (defined(Dades))

# declare NumElements = NumElements + 1;
#read (Dades, Alfai, Ai, Di, Fi)

#declare Taula[0][NumElements-1] = Alfai;
#declare Taula[1][NumElements-1] = Ai;
#declare Taula[2][NumElements  ] = Di;
#declare Taula[3][NumElements  ] = Fi;

#end 

/*
#fopen Calculs "calculs" write
#write (Calculs, X, "," , Y, "\n")
*/

/*}}}*/



/*{{{ Definicio de l'element 0 */

#declare Gruix=0.08;

#declare Element0= 

   union{
         cylinder{  <  0, -1, 0>, <0, 0, 0>, Gruix texture { Copper_Metal } }
        }
        
/*}}}*/
/*{{{ Definicio de l'element 1 */

#declare Element1= 

   union{
         cylinder{  <  0, 0, 0>, <2.0, 0, 0>, Gruix texture { Copper_Metal } }

         cylinder{  <0, 0, -0.5>,   <0, 0, 0.5>, 1.2*Gruix texture { Copper_Metal } }
        }
        
/*}}}*/
/*{{{ Definicio d'un element del mig */

#declare Element_del_mig= 

   union{
         cylinder{  <  0, 0, 0>, <2.0, 0, 0>, Gruix texture { Copper_Metal } }

         cylinder{<0, 0, -Gruix>, <0, 0, Gruix>, 1.2*Gruix  texture { Copper_Metal } }
        }
        
/*}}}*/
/*{{{ Definicio de l'extrem  */

#declare Extrem  = 

   union{

         cylinder{ <  0, 0, 0>, <1.5, 0, 0>, Gruix  }
         cone    { < 1.5, 0, 0 >   , Gruix    // Center and radius of one end
                  <  1.5+0.75*Gruix,0,0> , 0  // Center and radius of other end
                 }                
         cylinder{<0, 0, -Gruix>, <0, 0, Gruix>, 1.2*Gruix }
         texture { Copper_Metal }

        }
        
/*}}}*/


/*{{{ Escena */
// Escena 

// Base:
object {SRef}

object { Element0 }

// Sistema 1:

#declare Count = 1;

object { SRef 

           rotate    Taula[3][Count]*z
           translate Taula[2][Count]*z
           translate Taula[1][Count-1]*x
           rotate    Taula[0][Count-1]*x


       }

object { Element1

           rotate    Taula[3][Count]*z
           translate Taula[2][Count]*z
           translate Taula[1][Count-1]*x
           rotate    Taula[0][Count-1]*x


       }

/*{{{ Sistemes del mig:*/

#declare NumElement = NumElements-1;
#while (NumElement > 1)


union {
  object { Element_del_mig }
  object { SRef }

           #declare Count=NumElement;
	   #while (Count > 0)
           rotate    Taula[3][Count]*z
           translate Taula[2][Count]*z
           translate Taula[1][Count-1]*x
           rotate    Taula[0][Count-1]*x
	   #declare Count=Count-1;
	   #end
			 

       }

#declare NumElement=NumElement-1;

#end


/*}}}*/
       
/*{{{ Element extrem: */
// Sistema 3:

object { SRef  

           #declare Count=NumElements;
	   #while (Count > 0)
           rotate    Taula[3][Count]*z
           translate Taula[2][Count]*z
           translate Taula[1][Count-1]*x
           rotate    Taula[0][Count-1]*x
	   #declare Count=Count-1;
	   #end
			 
       }


object { Extrem

           #declare Count=NumElements;
	   #while (Count > 0)
           rotate    Taula[3][Count]*z
           translate Taula[2][Count]*z
           translate Taula[1][Count-1]*x
           rotate    Taula[0][Count-1]*x
	   #declare Count=Count-1;
	   #end
			 

       }

/*}}}*/

#declare Px =-1.2011115;
#declare Py = 0.3-0.5;
#declare Pz = 0.0;

  sphere { <Px, Py, Pz>, 0.3 texture { Copper_Metal }

         }



/*}}}*/

// END

