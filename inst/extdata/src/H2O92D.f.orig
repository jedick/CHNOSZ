*** H2O92 - Computes state, thermodynamic, transport, and electroststic 
***         properties of fluid H2O at T,[P,D] using equations and data
***         given by Haar et al. (1984), Levelt Sengers et al. (1983),
***         Johnson and Norton (1991), Watson et al. (1980), Sengers and
***         Kamgar-Parsi (1984), Sengers et al. (1984), Helgeson and Kirkham
***         (1974), Uematsu and Franck (1980), and Pitzer (1983). 
***
***********************************************************************
***
*** Author:     James W. Johnson
***             Earth Sciences Dept., L-219
***             Lawrence Livermore National Laboratory
***             Livermore, CA 94550
***             johnson@s05.es.llnl.gov
***
*** Abandoned:  8 November 1991
***
***********************************************************************
*
*   specs  - Input unit, triple point, saturation, and option specs:
*
*****        it, id, ip, ih, itripl, isat, iopt, useLVS, epseqn, icrit;
*
*            note that the returned value of isat may differ from
*            its input value and that icrit need not be specified
*            prior to invocation.
*
*
*   states - State variables: 
*
*****          temp, pres, dens(1), dens(2);
*
*            note that the first three of these must be specified prior
*            to invocation and that, in the case of saturation, vapor
*            density is returned in dens(1), liquid in dens(2).
*
*
*   props  - Thermodynamic, transport, electrostatic, and combined 
*            property values:
*
*****        A, G, S, U, H, Cv, Cp, Speed, alpha, beta, diel, visc,
*****        tcond, surten, tdiff, Prndtl, visck, albe, 
*****        ZBorn, YBorn, QBorn, daldT, XBorn
*
*
*   error  - LOGICAL argument that indicates success ("FALSE") or
*            failure ("TRUE") of the call, the latter value in 
*            response to out-of-bounds specs or states variables.
*
***********************************************************************

      SUBROUTINE H2O92(specs,states,props,error)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP  = 23, NPROP2 = 46)

      INTEGER  specs(10)
      DOUBLE PRECISION  states(4), props(NPROP2), Dens(2),
     1                  wpliq(NPROP), wprops(NPROP)
      LOGICAL           crtreg, valid, error, useLVS

      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
      COMMON /wpvals/ wprops, wpliq

      SAVE


      CALL unit(specs(1),specs(2),specs(3),specs(4),specs(5))

      IF (.NOT. (valid(specs(1),specs(2),specs(3),specs(4),specs(5),
     1                 specs(6),specs(7),specs(8),specs(9),
     2                 states(1),states(2),states(3)))) THEN
           error = .TRUE.
           RETURN
      ELSE
           error = .FALSE.
      END IF

      IF (crtreg(specs(6),specs(7),specs(1),
     1           states(1),states(2),states(3))) THEN
           specs(10) = 1
           useLVS = (specs(8) .EQ. 1)
      ELSE
           specs(10) = 0
           useLVS = .FALSE.
      END IF


      IF (useLVS) THEN
           Dens(1) = states(3)
           CALL LVSeqn(specs(6),specs(7),specs(5),
     1                 states(1),states(2),Dens,specs(9))
           Dens(1) = Dens(1) / 1.0d3
           IF (specs(6) .EQ. 1) THEN
                Dens(2) = Dens(2) / 1.0d3
           END IF
      ELSE
           Dens(1) = states(3) / 1.0d3
           CALL HGKeqn(specs(6),specs(7),specs(5),
     1                 states(1),states(2),Dens,specs(9))
      END IF

      CALL load(1,wprops,props)

      IF (specs(6) .EQ. 1) THEN
           tempy = Dens(1)
           Dens(1) = Dens(2)
           Dens(2) = tempy
           CALL load(2,wpliq,props)
      END IF
      
      states(1) = TdegUS(specs(1),states(1))
      states(2) = states(2) * fp
      states(3) = Dens(1) / fd

      IF (specs(6) .EQ. 1) THEN
           states(4) = Dens(2) / fd
      END IF

      RETURN
      END

************************************************************************

*** valid - Returns "TRUE" if unit and equation specifications 
*           are valid and input state conditions fall within
*           the HGK equation's region of validity;
*           returns "FALSE" otherwise.

      LOGICAL FUNCTION valid(it,id,ip,ih,itripl,isat,iopt,
     1                       useLVS,epseqn,Temp,Pres,Dens)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER  useLVS, epseqn
      LOGICAL  valspc, valTD, valTP

      COMMON /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
      COMMON /crits/  Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /tpoint/ Utr, Str, Htr, Atr, Gtr,
     1                Ttr, Ptripl, Dltrip, Dvtrip
      COMMON /HGKbnd/ Ttop, Tbtm, Ptop, Pbtm, Dtop, Dbtm
      COMMON /liqice/ sDli1, sPli1, sDli37, sPli37, sDIB30,
     1                Tli13, Pli13, Dli13, TnIB30, DnIB30

      SAVE


*** ensure validity of input specifications
      IF (.NOT. valspc(it,id,ip,ih,itripl,isat,iopt,
     1                 useLVS,epseqn)) THEN
           valid = .FALSE.
           RETURN
      END IF

*** convert to  degC, bars, g/cm3 ***
      T = TdegK(it,Temp) - 273.15d0
      D = Dens * fd
      P = Pres / fp * 1.0d1
      Ttripl = Ttr - 273.15d0
      Tcrit = Tc - 273.15d0
      Pcrit = Pc * 1.0d1

      IF (isat .EQ. 0) THEN
           IF (iopt .EQ. 1) THEN
                valid = valTD(T,D,isat,epseqn)
           ELSE
                valid = valTP(T,P)
           END IF
      ELSE
           IF (iopt .EQ. 1) THEN
                valid = ((T+FPTOL .GE. Ttripl) .AND. 
     1                   (T-FPTOL .LE. Tcrit))
           ELSE
                valid = ((P+FPTOL .GE. Ptripl) .AND. 
     1                   (P-FPTOL .LE. Pcrit))
           END IF
      END IF

      RETURN
      END

*****************************************************************

*** valspc - Returns "TRUE" if  it, id, ip, ih, itripl, isat, iopt,
*            useLVS, and epseqn values all define valid input;
*            returns "FALSE" otherwise.

      LOGICAL FUNCTION valspc(it,id,ip,ih,itripl,isat,iopt,
     1                        useLVS,epseqn)

      INTEGER  useLVS, epseqn

      SAVE

      
      valspc = (1 .LE. it)     .AND. (it     .LE. 4) .AND.
     1         (1 .LE. id)     .AND. (id     .LE. 4) .AND.
     2         (1 .LE. ip)     .AND. (ip     .LE. 5) .AND.
     3         (1 .LE. ih)     .AND. (ih     .LE. 6) .AND.
     4         (0 .LE. itripl) .AND. (itripl .LE. 1) .AND.
     5         (0 .LE. isat)   .AND. (isat   .LE. 1) .AND.
     6         (1 .LE. iopt)   .AND. (iopt   .LE. 2) .AND.
     7         (0 .LE. useLVS) .AND. (useLVS .LE. 1) .AND. 
     8         (1 .LE. epseqn) .AND. (epseqn .LE. 5)

      RETURN
      END

*****************************************************************

*** valTD - Returns "TRUE" if  T-D  defines liquid or vapor H2O
*           within validity limits of the HGK equation of state;
*           returns "FALSE" otherwise.

      LOGICAL FUNCTION valTD(T,D,isat,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER epseqn

      COMMON /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
      COMMON /RTcurr/ rt
      COMMON /crits/  Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /tpoint/ Utr, Str, Htr, Atr, Gtr,
     1                Ttr, Ptripl, Dltrip, Dvtrip
      COMMON /HGKbnd/ Ttop, Tbtm, Ptop, Pbtm, Dtop, Dbtm
      COMMON /liqice/ sDli1, sPli1, sDli37, sPli37, sDIB30,
     1                Tli13, Pli13, Dli13, TnIB30, DnIB30
      COMMON /coefs/  a(20), q(20), x(11)
      COMMON /satur/  Dliq, Dvap, DH2O, iphase

      SAVE
   
      EQUIVALENCE  (TmnLVS, x(1))


      IF ((T-FPTOL .GT. Ttop) .OR. (T+FPTOL .LT. Tbtm) .OR. 
     1    (D-FPTOL .GT. Dtop) .OR. (D+FPTOL .LT. Dbtm)) THEN
           valTD = .FALSE.
           RETURN
      END IF

      Tcrit = Tc - 273.15d0
      Ttripl = Ttr - 273.15d0

      IF ((T+FPTOL .GE. Tcrit)  .OR.
     1   ((T .GE. TnIB30) .AND. (D .GE. Dltrip))) THEN
           Dlimit = sDIB30 * (T-TnIB30) + Dtop
           valTD  = (D-FPTOL .LE. Dlimit)
      ELSE
           IF (D-FPTOL .LE. Dltrip) THEN
                IF (T .GE. Ttripl) THEN
                     valTD = .TRUE.
                     Tk = T + 273.15d0
                     IF (Tk .LT. TmnLVS) THEN
                          rt = gascon * Tk
                          CALL pcorr(0,Tk,Ps,Dl,Dv,epseqn)
                     ELSE
                          istemp = 1
                          DH2O = 0.0d0
                          P = Pfind(istemp,Tk,DH2O)
                          CALL denLVS(istemp,Tk,P)
                          Dv = Dvap / 1.0d3
                          Dl = Dliq / 1.0d3
                     END IF
                     IF ((D .GE. Dv) .AND. (D. LE. Dl)) THEN
                          isat = 1
                     END IF
                ELSE
                     P = Psublm(T)
                     PMPa = P / 1.0d1
                     Tk = T + 273.15d0
                     Dguess = PMPa / Tk / 0.4d0
                     rt = gascon * Tk
                     CALL bb(Tk)
                     CALL denHGK(Dsublm,PMPa,Dguess,Tk,dPdD)
                     valTD = (D-FPTOL .LE. Dsublm)
                END IF
           ELSE
                IF (D .LE. Dli13) THEN
                     Dlimit = sDli1 * (T-Tli13) + Dli13
                     valTD = (D+FPTOL .GE. Dlimit)
                ELSE
                     Dlimit = sDli37 * (T-Tli13) + Dli13
                     valTD = (D-FPTOL .LE. Dlimit)
                END IF
           END IF
      END IF

      RETURN
      END

*****************************************************************

*** valTP - Returns "TRUE" if  T-P  defines liquid or vapor H2O
*           within validity limits of the HGK equation of state;
*           returns "FALSE" otherwise.

      LOGICAL FUNCTION valTP(T,P)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      COMMON /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
      COMMON /crits/  Tcrit, rhoC, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /tpoint/ Utr, Str, Htr, Atr, Gtr,
     1                Ttr, Ptripl, Dltrip, Dvtrip
      COMMON /HGKbnd/ Ttop, Tbtm, Ptop, Pbtm, Dtop, Dbtm
      COMMON /liqice/ sDli1, sPli1, sDli37, sPli37, sDIB30,
     1                Tli13, Pli13, Dli13, TnIB30, DnIB30

      SAVE


      IF ((T-FPTOL .GT. Ttop) .OR. (T+FPTOL .LT. Tbtm) .OR.
     1    (P-FPTOL .GT. Ptop) .OR. (P+FPTOL .LT. Pbtm)) THEN
           valTP = .FALSE.
           RETURN
      ELSE
           valTP = .TRUE.
      END IF

      IF (P .GE. Pli13) THEN
           Plimit = sPli37 * (T-Tli13) + Pli13
           valTP = (P-FPTOL .LE. Plimit)
      ELSE
           IF (P .GE. Ptripl) THEN
                Plimit = sPli1 * (T-Tli13) + Pli13
                valTP = (P+FPTOL .GE. Plimit)
           ELSE
                Psubl = Psublm(T)
                valTP = (P-FPTOL .LE. Psubl)
           END IF
      END IF

      RETURN
      END

*****************************************************************

*** Psublm - Returns  Psublimation(T)  computed from the 
*            equation given by  Washburn (1924): Monthly
*            Weather Rev., v.52, pp.488-490.

      DOUBLE PRECISION FUNCTION Psublm(Temp)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      SAVE

      
      T = Temp + 2.731d2

      PmmHg = power(1.0d1, (-2.4455646d3/T + 8.2312d0*DLOG10(T) - 
     1              1.677006d-2*T + 1.20514d-5*T*T - 6.757169d0))

*** convert mmHg to bars ***
      Psublm = PmmHg * 1.33322d-3

      RETURN
      END

************************************************************************

*** HGKcon - Constant parameters for the H2O equation of state 
*            given by  Haar, Gallagher, & Kell (1984): 
*            bp, bq     = b(j), B(j) from Table A.1, p.272
*            g1, g2, gf = alpha, beta, gamma from eq (A.2), p.272
*            g, ii, jj  = g(i), k(i), l(i) from eq (A.5), p.272.
*            Note that  tz < tcHGK.
*                 Tolerence limits required in various real & inexact
*            comparisons are set and stored in COMMON /tolers/.


      BLOCK DATA HGKcon

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
      COMMON /nconst/ g(40), ii(40), jj(40), nc
      COMMON /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
      COMMON /bconst/ bp(10), bq(10)
      COMMON /addcon/ atz(4), adz(4), aat(4), aad(4)
      COMMON /HGKcrt/ tcHGK, dcHGK, pcHGK
      COMMON /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
      COMMON /HGKbnd/ Ttop, Tbtm, Ptop, Pbtm, Dtop, Dbtm
      COMMON /liqice/ sDli1, sPli1, sDli37, sPli37, sDIB30,
     1                Tli13, Pli13, Dli13, TnIB30, DnIB30
      COMMON /tpoint/ Utripl, Stripl, Htripl, Atripl, Gtripl, 
     1                Ttripl, Ptripl, Dltrip, Dvtrip

      SAVE

      DATA    Ttripl, Ptripl, Dltrip, Dvtrip
     1     /  2.7316d2,
     2        0.611731677193563186622762580414d-2,
     3        0.999778211030936587977889295063d0,
     4        0.485467583448287303988319166423d-5 /

      DATA   Ttop,    Tbtm,    Ptop,    Pbtm,    Dtop,    Dbtm
     1    /  2.25d3, -2.0d1,   3.0d4,   1.0d-3, 
     2       0.138074666423686955066817336896d1,
     3       0.858745555396173972667420987465d-7 /

      DATA  sDli1, sPli1, sDli37, sPli37, sDIB30,
     1      Tli13, Pli13, Dli13, TnIB30, DnIB30
     2    / -0.584797401732178547634910059828d-2,
     3      -0.138180804975562958027981345769d3,
     4       0.183244000000000000000000000007d-2,
     5       0.174536874999999999999999999995d3,
     6      -0.168375439429928741092636579574d-3,
     7      -0.15d2,
     8       0.20741d4,
     9       0.108755631570602617113573577945d1,
     1       0.145d3,
     1       0.102631640581853166397515716306d1 /


      DATA    TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL 
     1     / 1.0d-6, 1.0d-6, 1.0d-9, 1.0d-5, -673.5d0, 1.0d-7 /

      DATA tcHGK, dcHGK, pcHGK / .647126d3, .322d3, .22055d2 /

      DATA atz /.64d3,  .64d3,  .6416d3, .27d3/
      DATA adz /.319d0, .319d0, .319d0,  .155d1/
      DATA aat /.2d5,   .2d5,   .4d5,    .25d2/
      DATA aad /.34d2,  .4d2,   .3d2,    .105d4/

      DATA wm, gascon, tz, aa, uref, sref
     1    /  .1801520000d2,  .46152200d0,  .647073d3,  .1d1, 
     2      -.4328455039d4,  .76180802d1 /

      DATA g1, g2, gf /.11d2, .44333333333333d2, .35d1/

      DATA bp / .7478629d0,  -.3540782d0,  2*.0d0,  .7159876d-2, 
     1          .0d0,        -.3528426d-2, 3*.0d0/

      DATA bq / .11278334d1,  .0d0, -.5944001d0, -.5010996d1, .0d0,
     1          .63684256d0,  4*.0d0/

      DATA nc / 36 /

      DATA g /-.53062968529023d3,  .22744901424408d4,  .78779333020687d3
     1,       -.69830527374994d2,  .17863832875422d5, -.39514731563338d5
     2,        .33803884280753d5, -.13855050202703d5, -.25637436613260d6
     3,        .48212575981415d6, -.34183016969660d6,  .12223156417448d6
     4,        .11797433655832d7, -.21734810110373d7,  .10829952168620d7
     5,       -.25441998064049d6, -.31377774947767d7,  .52911910757704d7
     6,       -.13802577177877d7, -.25109914369001d6,  .46561826115608d7
     7,       -.72752773275387d7,  .41774246148294d6,  .14016358244614d7
     8,       -.31555231392127d7,  .47929666384584d7,  .40912664781209d6
     9,       -.13626369388386d7,  .69625220862664d6, -.10834900096447d7
     a,       -.22722827401688d6,  .38365486000660d6,  .68833257944332d4
     b,        .21757245522644d5, -.26627944829770d4, -.70730418082074d5
     c,       -.22500000000000d0, -.16800000000000d1
     d,        .5500000000000d-1, -.93000000000000d2/

      DATA ii / 4*0, 4*1, 4*2, 4*3, 4*4, 4*5, 4*6, 4*8, 2*2, 0, 4, 
     1          3*2, 4/

      DATA jj / 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 
     1          2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7,
     2          1, 4, 4, 4, 0, 2, 0, 0/

      END

*********************************************************************

*** LVScon - Constant parameters for the H2O critical region equation 
*            of state given by  Levelt Sengers, Kamgar-Parsi, Balfour,
*            & Sengers (1983).

      BLOCK DATA LVScon

      IMPLICIT DOUBLE PRECISION  (a-h,o-z)

      COMMON /crits/ Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /coefs/ a(20), q(20), x(11)

      SAVE
     
      DATA   Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon
     1       /   647.067d0, 322.778d0, 22.046d0, 
     2           0.034070660379837018423130834983d0, 22046.0d0,
     3           0.034070660379837018423130834983d3,
     4           0.000000327018783663660700780197d0 /

      DATA a /  -0.017762d0,  5.238000d0,  0.000000d0, -2.549150d1,  
     1           6.844500d0,  0.325000d0,  1.440300d0,  0.000000d0, 
     2           1.375700d0,  2.366660d1,  4.820000d0,  0.294200d0,
     3          -1.123260d1, -2.265470d1, -1.788760d1, -4.933200d0,
     4           1.109430391161019373812391218008d0,
     5          -1.981395981400671095301629432211d0,
     6           0.246912528778663959151808173743d0,
     7          -0.843411332867484343974055795059d0 / 

      DATA q /  -0.006000d0, -0.003000d0,  0.000000d0,  6.470670d2,
     1           3.227780d2,  2.204600d1,  0.267000d0, -1.600000d0,
     2           0.491775937675717720291497417773d0,    0.108500d0,
     3           0.586534703230779473334597524774d0,
     4          -1.026243389120214352553706598564d0,
     5           0.612903225806451612903225804745d0,    0.500000d0,
     6          -0.391500d0,  0.825000d0,  0.741500d0,
     7           0.103245882826119154987166286332d0,
     8           0.160322434159191991394857495360d0,
     9          -0.169859514687100893997445721324d0 /

      DATA x /   6.430000d2,  6.453000d2,  6.950000d2,
     1           1.997750d2,  4.200400d2,
     2           2.09945691135940719075293945960d1,
     3           2.15814057875264119875397458907d1,
     4           3.0135d1, 4.0484d1,
     5           .175777517046267847932127026995d0,
     6           .380293646126229135059562456934d0 /


*     EQUIVALENCE (cc,     a(1) ),  (pointA, q(1) ),  (Tmin1,  x(1)),
*    1            (p3,     a(2) ),  (pointB, q(2) ),  (Tmin2,  x(2)),
*    2            (delroc, a(3) ),  (delpc,  q(3) ),  (Tmax,   x(3)),
*    3            (p2,     a(4) ),  (Tc,     q(4) ),  (Dmin,   x(4)),
*    4            (p1,     a(5) ),  (rhoc,   q(5) ),  (Dmax,   x(5)),
*    5            (beta,   a(6) ),  (Pc,     q(6) ),  (Pmin1,  x(6)),
*    6            (xko,    a(7) ),  (dPcdTc, q(7) ),  (Pmin2,  x(7)),
*    7            (delTc,  a(8) ),  (slopdi, q(8) ),  (Pmax1,  x(8)),
*    8            (besq,   a(9) ),  (p11,    q(9) ),  (Pmax2,  x(9)),
*    9            (aa,     a(10)),  (alpha,  q(10)),  (sl1,    x(10)),
*    0            (delta,  a(11)),  (p00,    q(11)),  (sl2,    x(11)),
*    1            (k1,     a(12)),  (p20,    q(12)),
*    2            (muc,    a(13)),  (p40,    q(13)),
*    3            (mu1,    a(14)),  (deli,   q(14)),
*    4            (mu2,    a(15)),  (alh1,   q(15)),
*    5            (mu3,    a(16)),  (beti,   q(16)),
*    6            (s00,    a(17)),  (gami,   q(17)),
*    7            (s20,    a(18)),  (p01,    q(18)),
*    8            (s01,    a(19)),  (p21,    q(19)),
*    9            (s21,    a(20)),  (p41,    q(20))

      END

*******************************************************************

*** unit - Sets internal parameters according to user-specified 
*          choice of units.  Internal program units are degK(T),
*          and gm/cm**3(D); all other properties are computed in
*          dimensionless form and dimensioned at output time.
*          NOTE:  conversion factors for j/g ---> cal/(g,mole)
*          (ffh (4 & 5)) are consistent with those given in 
*          Table 1, Helgeson & Kirkham (1974a) for thermal calories, 
*          and differ slightly with those given by Haar et al (1984) 
*          for international calories.  

      SUBROUTINE unit(it,id,ip,ih,itripl)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION  fft(4), ffd(4), ffvd(4), ffvk(4), 
     1                  ffs(4), ffp(5), ffh(6),
     2                  ffst(4), ffcd(4), ffch(6)

      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc

      SAVE

      DATA fft  /1.0d0,  1.0d0, 0.555555556d0, 0.555555556d0 /
      DATA ffd  /1.0d-3, 1.0d0, 1.80152d-2,    1.6018d-2/
      DATA ffvd /1.0d0,  1.0d1, 0.555086816d0, 0.671968969d0 /
      DATA ffvk /1.0d0,  1.0d4, 1.0d4,         1.076391042d1 /
      DATA ffs  /1.0d0,  1.0d2, 1.0d2,         3.280833d0 /
      DATA ffp  /1.0d0,  1.0d1, 9.869232667d0, 1.45038d2,   1.01971d1/
      DATA ffh  /1.0d0,  1.0d0, 1.80152d1,     2.3901d-1,
     1           4.305816d0, 4.299226d-1/
      DATA ffst /1.0d0,  1.0d3,  0.555086816d2, 0.2205061d1 /
      DATA ffcd /1.0d0,  1.0d-2, 1.0d-2,        0.3048d0 /
      DATA ffch /1.0d-3, 1.0d0,  1.0d0,         0.23901d0,
     1           0.23901d0, 0.947244d-3 /


      ft  = fft(it)
      fd  = ffd(id)
      fvd = ffvd(id)
      fvk = ffvk(id)
      fs  = ffs(id)
      fp  = ffp(ip)
      fh  = ffh(ih)
      fst = ffst(id)
      fc  = ffcd(id) * ffch(ih)

      IF (itripl .EQ. 1)  CALL tpset

      RETURN

      END

***********************************************************************

*** crtreg - Returns "TRUE" if input state conditions fall within 
*            the critical region of H2O; otherwise returns "FALSE".
*            T, P, D, input in user-specified units, are returned in
*            degK, MPa, kg/m3.

      LOGICAL FUNCTION crtreg(isat,iopt,it,T,P,D)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      LOGICAL  llim, ulim

      COMMON /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /coefs/ a(20), q(20), x(11)
      COMMON /units/ ft, fd, fvd, fvk, fs, fp, fh, fst, fc

      SAVE

      EQUIVALENCE (Tmin1,  x(1)),  (Tmin2, x(2)),  (Tmax,  x(3)),
     1            (Dmin,   x(4)),  (Dmax,  x(5)),
     2            (Pbase1, x(6)),  (Pbase2,x(7)),
     3            (PTmins, x(10)), (PTmaxs,x(11))


      T = TdegK(it,T)
      IF (isat .EQ. 0) THEN
           IF (iopt .EQ. 1) THEN
                D = D * fd * 1.0d3
                crtreg = ((T .GE. Tmin1) .AND. (T .LE. Tmax) .AND.
     1                    (D .GE. Dmin) .AND. (D .LE. Dmax))
           ELSE
                P = P / fp
                IF ((T .LT. Tmin1) .OR. (T .GT. Tmax)) THEN
                     crtreg = .FALSE.
                ELSE
                     Pmin = Pbase1 + PTmins * (T - Tmin1)
                     Pmax = Pbase2 + PTmaxs * (T - Tmin2)
                     llim = (P .GE. Pmin)
                     ulim = (P .LE. Pmax)
                     IF (llim .AND. ulim) THEN
                          crtreg = .TRUE.
                     ELSE
                          IF (llim .AND. (T .LE. Tmin2)) THEN
                               isat1 = 1
                               ddummy = 0.0d0
                               Pstest = Pfind(isat1,T,ddummy)
                               crtreg = (P .LE. Pstest)
                          ELSE
                               crtreg = .FALSE.
                          END IF
                     END IF
                END IF
           END IF
      ELSE
           IF (iopt .EQ. 1) THEN
                crtreg = (T .GE. Tmin1)
           ELSE
                P = P / fp
                crtreg = (P .GE. Pbase1)
           END IF
      END IF

      RETURN
      END

*********************************************************************

*** HGKeqn - Computes thermodynamic and transport properties of 
*            of H2O from the equation of state given by  
*            Haar, Gallagher, & Kell (1984).

      SUBROUTINE HGKeqn(isat,iopt,itripl,Temp,Pres,Dens,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP = 23)

      INTEGER epseqn
      DOUBLE PRECISION  Dens(2), wprops(NPROP), wpliq(NPROP)

      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
      COMMON /wpvals/ wprops, wpliq
      COMMON /RTcurr/ rt

      SAVE

      rt = gascon * Temp

      CALL HGKsat(isat,iopt,itripl,Temp,Pres,Dens,epseqn)

      IF (isat .EQ. 0) THEN
           CALL bb(Temp)
           CALL calcv3(iopt,itripl,Temp,Pres,Dens(1),epseqn)
           CALL thmHGK(Dens(1),Temp)
           CALL dimHGK(isat,itripl,Temp,Pres,Dens(1),epseqn)
      ELSE
           DO 10  i=1,NPROP
 10             wpliq(i) = wprops(i)
           CALL dimHGK(2,itripl,Temp,Pres,Dens(2),epseqn)
      END IF

      RETURN
      END

*****************************************************************

*** HGKsat - If  isat=1, computes  Psat(T) or Tsat(P) (iopt=1,2),
*            liquid and vapor densities, and associated
*            thermodynamic and transport properties.
*            If  isat=0, checks whether  T-D or T-P (iopt=1,2)
*            falls on or within  TOL  of the liquid-vapor 
*            surface; if so, sets isat <- 1 and computes 
*            properties.
   
      SUBROUTINE HGKsat(isat,iopt,itripl,Temp,Pres,Dens,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      
      DOUBLE PRECISION  Dens(2)
      INTEGER  epseqn

      COMMON /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
      COMMON /HGKcrt/ tcHGK, dcHGK, pcHGK
      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
      COMMON /tpoint/ Utr, Str, Htr, Atr, Gtr, 
     1                Ttripl, Ptripl, Dltrip, Dvtrip
      COMMON /crits/  Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon

      SAVE


      IF (isat .EQ. 1) THEN
           IF (iopt .EQ. 1) THEN
                CALL pcorr(itripl,Temp,Pres,Dens(1),Dens(2),epseqn)
           ELSE
                CALL tcorr(itripl,Temp,Pres,Dens(1),Dens(2),epseqn)
           END IF
      ELSE
           IF ((Temp .GT. Tc) .OR. (Temp .LT. Ttripl) .OR.
     1        ((iopt .EQ. 2) .AND. (Pres .GT. Pc))) THEN
                RETURN
           ELSE
                CALL pcorr(itripl,Temp,Ptemp,dltemp,dvtemp,epseqn)
                IF (((iopt .EQ. 2) .AND. 
     1              (DABS(Pres-Ptemp) .LE. PTOL)) .OR.
     2              ((iopt .EQ. 1) .AND. 
     3              ((DABS(Dens(1)-dltemp) .LE. DTOL) .OR.
     4              (DABS(Dens(1)-dvtemp) .LE. DTOL)))) THEN
                          isat = 1
                          Pres = Ptemp
                          Dens(1) = dltemp
                          Dens(2) = dvtemp
                END IF
           END IF
      END IF

      RETURN
      END

************************************************************************

*** calcv3 - Compute the dependent state variable.

      SUBROUTINE calcv3(iopt,itripl,Temp,Pres,Dens,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER  epseqn

      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
      COMMON /qqqq/   q0, q5
      COMMON /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
      COMMON /fcts/   ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, dpdd,
     1                cjtt, cjth
      COMMON /RTcurr/ rt

      SAVE


      IF (iopt .EQ. 1) THEN
           CALL resid(Temp,Dens)
           CALL base(Dens,Temp)
           CALL ideal(Temp)
           Pres  = rt * Dens * z + q0
      ELSE
           IF (Temp .LT. tz) THEN  
                CALL pcorr(itripl,Temp,ps,dll,dvv,epseqn)
           ELSE
                ps   = 2.0d4
                dll  = 0.0d0
           END IF
           IF (Pres .GT. ps) THEN
                dguess = dll
           ELSE
                dguess = Pres / Temp / 0.4d0
           END IF

           CALL denHGK(Dens,Pres,dguess,Temp,dpdd)
           CALL ideal(Temp)
      END IF

      RETURN
      END

******************************************************************************

*** thmHGK - Computes thermodynamic functions in dimensionless
*            units from the HGK equation of state:  Helmholtz, Gibbs, 
*            internal energy, and enthalpy functions (ad, gd, ud, hd) are 
*            per RT; entropy and heat capacities (sd, cvd, cpd) are per R.

      SUBROUTINE thmHGK(d,t)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, y, uref, sref
      COMMON /qqqq/   qp, qdp
      COMMON /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
      COMMON /resf/   ar, gr, sr, ur, hr, cvr, dpdtr
      COMMON /idf/    ai, gi, si, ui, hi, cvi, cpi
      COMMON /fcts/   ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, dpdd,
     1                cjtt, cjth
      COMMON /RTcurr/ rt

      SAVE


      z    = zb + qp/rt/d
      dpdd = rt * (zb + y * dzb) + qdp
      ad   = ab + ar + ai - uref/t + sref
      gd   = ad + z
      ud   = ub + ur + ui - uref/t
      dpdt = rt * d * dpdtb + dpdtr
      cvd  = cvb + cvr + cvi
      cpd  = cvd + t*dpdt*dpdt/(d*d*dpdd*gascon)
      hd   = ud + z
      sd   = sb + sr + si - sref
      dvdt = dpdt / dpdd / d / d
      cjtt = 1.0d0 / d - t * dvdt
      cjth = -cjtt / cpd / gascon

      RETURN

      END

*************************************************************************

*** bb - Computes molecular parameters b, the "excluded volume" 
*        (eq A.3), and B, the second virial coefficient (eq A.4),
*        in cm3/g (b1,b2) and their first and second derivatives 
*        with respect to temperature (b1t,b1tt,b2t,b2tt).

      SUBROUTINE bb(t)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION v(10)

      COMMON /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
      COMMON /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
      COMMON /bconst/ bp(10), bq(10)

      SAVE


      v(1) = 1.0d0
     
      DO 2 i=2,10
 2         v(i) = v(i-1) * tz / t

      b1   = bp(1) + bp(2) * DLOG(1.0 / v(2))
      b2   = bq(1)
      b1t  = bp(2) * v(2) / tz
      b2t  = 0.0d0
      b1tt = 0.0d0
      b2tt = 0.0d0

      DO 4 i=3,10
           b1   = b1   + bp(i) * v(i-1)
           b2   = b2   + bq(i) * v(i-1)
           b1t  = b1t  - (i-2) * bp(i) * v(i-1) / t
           b2t  = b2t  - (i-2) * bq(i) * v(i-1) / t
           b1tt = b1tt + bp(i) * (i-2)*(i-2) * v(i-1) / t / t
 4         b2tt = b2tt + bq(i) * (i-2)*(i-2) * v(i-1) / t / t

      b1tt = b1tt - b1t / t
      b2tt = b2tt - b2t / t

      RETURN

      END

***********************************************************************

*** base - Computes Abase, Gbase, Sbase, Ubase, Hbase, Cvbase
*          -- all per RT (dimensionless) --  as well as Pbase & dP/dT
*          -- both per (DRT) -- for the base function (ab, gb, sb, ub,
*          hb, cvb, pb, dpdtb).  See Haar, Gallagher & Kell (1979), eq(1).

                
      SUBROUTINE base(d,t)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      
      COMMON /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
      COMMON /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
      COMMON /aconst/ wm, gascon, tz, a, z, dz, y, uref, sref

      SAVE


      y     = .25d0 * b1 * d
      x     = 1.0d0 - y
      z0    = (1.0d0 + g1*y + g2*y*y) / (x*x*x)
      z     = z0 + 4.0d0*y*(b2/b1 - gf)
      dz0   = (g1 + 2.0d0*g2*y)/(x*x*x) + 
     1        3.0d0*(1.0d0 + g1*y + g2*y*y)/(x*x*x*x)
      dz    = dz0 + 4.0d0*(b2/b1 - gf)

      pb    = z

      ab    = -DLOG(x) - (g2 - 1.0d0)/x + 28.16666667d0/x/x +
     1         4.0d0*y*(b2/b1 - gf) + 15.166666667d0 +
     2         DLOG(d*t*gascon/.101325d0)
      gb    = ab + z
      ub    = -t*b1t*(z - 1.0d0 - d*b2)/b1 - d*t*b2t
      sb    = ub - ab
      hb    = z + ub

      bb2tt = t * t * b2tt
      cvb   = 2.0d0*ub + (z0 - 1.0d0)*(((t*b1t/b1)*(t*b1t/b1)) - 
     1        t*t*b1tt/b1) - d*(bb2tt - gf*b1tt*t*t) - 
     2        (t*b1t/b1)*(t*b1t/b1)*y*dz0
      
      dpdtb = pb/t + d*(dz*b1t/4.0d0 + b2t - b2/b1*b1t)

      RETURN

      END

***********************************************************************

*** resid - Computes residual contributions to pressure (q), the 
*           Helmloltz function (ar) , dP/dD (q5), the Gibbs function 
*           (gr), entropy (sr), internal energy (ur), enthalpy (hr), 
*           isochoric heat capacity (cvr), and dP/dT.  The first 36 
*           terms of the residual function represent a global
*           least-squares fit to experimental data outside the 
*           critical region, terms 37-39 affect only the immediate
*           vicinity of the critical point, and the last term (40)
*           contributes only in the high pressure, low temperature
*           region.

      SUBROUTINE resid(t,d)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION qr(11), qt(10), qzr(9), qzt(9)

      COMMON /resf/   ar, gr, sr, ur, hr, cvr, dpdtr
      COMMON /qqqq/   q, q5
      COMMON /nconst/ g(40), ii(40), jj(40), n
      COMMON /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
      COMMON /addcon/ atz(4), adz(4), aat(4), aad(4)
      COMMON /RTcurr/ rt
      COMMON /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL

      SAVE

      EQUIVALENCE (qr(3), qzr(1)), (qt(2), qzt(1))


      qr(1) = 0.0d0
      q5    = 0.0d0
      q     = 0.0d0
      ar    = 0.0d0
      dadt  = 0.0d0
      cvr   = 0.0d0
      dpdtr = 0.0d0

      e     = DEXP(-aa * d)
      q10   = d * d * e
      q20   = 1.0d0 - e
      qr(2) = q10
      v     = tz / t
      qt(1) = t / tz

      DO 4 i=2,10
           qr(i+1) = qr(i) * q20
 4         qt(i)   = qt(i-1) * v

      DO 10 i=1,n
           k     = ii(i) + 1
           l     = jj(i)
           zz    = k
           IF (k .EQ. 1) THEN
                qp    = g(i) * aa * qr(2) * qzt(l)
           ELSE
                qp    = g(i) * aa * qzr(k-1) * qzt(l)
           END IF
           q     = q + qp
           q5    = q5 + aa*(2.0/d - aa*(1.0 - e*(k-1)/q20))*qp
           ar    = ar + g(i)*qzr(k)*qzt(l)/q10/zz/rt
           dfdt  = power(q20,DBLE(k))*(1-l)*qzt(l+1)/tz/k
           d2f   = l * dfdt
           dpt   = dfdt*q10*aa*k/q20
           dadt  = dadt  + g(i)*dfdt
           dpdtr = dpdtr + g(i)*dpt
 10        cvr   = cvr   + g(i)*d2f/gascon

      qp  = 0.0d0
      q2a = 0.0d0

      DO 20 j=37,40
           IF (g(j) .EQ. 0.0d0) GO TO 20
           k     = ii(j)
           km    = jj(j)
           ddz   = adz(j-36)
           del   = d/ddz - 1.0d0
           IF (DABS(del) .LT. 1.0d-10)  del = 1.0d-10
           ex1   = -aad(j-36) * power(del,DBLE(k))
           IF (ex1 .LT. EXPTOL) THEN
                dex = 0.0d0
           ELSE
                dex = DEXP(ex1)  * power(del,DBLE(km))
           END IF
           att   = aat(j-36)
           tx    = atz(j-36)
           tau   = t/tx - 1.0d0
           ex2   = -att * tau * tau
           IF (ex2 .LE. EXPTOL) THEN
                tex = 0.0d0
           ELSE
                tex = DEXP(ex2)
           END IF
           q10   = dex * tex
           qm    = km/del - k*aad(j-36)*power(del,DBLE(k-1))
           fct   = qm * d*d * q10 / ddz
           q5t   = fct*(2.0d0/d + qm/ddz) - (d/ddz)*(d/ddz)*q10 * 
     1             (km/del/del + k*(k-1)*aad(j-36) *
     2             power(del,DBLE(k-2)))
           q5    = q5 + q5t*g(j)
           qp    = qp + g(j)*fct
           dadt  = dadt  - 2.0d0*g(j)*att*tau* q10 /tx
           dpdtr = dpdtr - 2.0d0*g(j)*att*tau* fct /tx
           q2a   = q2a + t*g(j)*(4.0d0*att*ex2 + 2.0d0*att)*q10/tx/tx
           ar    = ar  + q10*g(j)/rt
 20        CONTINUE

      sr  = -dadt / gascon
      ur  = ar + sr
      cvr = cvr + q2a/gascon
      q   = q + qp

      RETURN

      END

************************************************************************

*** ideal - Computes thermodynamic properties for H2O in the 
*           ideal gas state using equations given by Woolley (1979).

      SUBROUTINE ideal(t)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION c(18)

      COMMON /idf/ ai, gi, si, ui, hi, cvi, cpi

      SAVE

      DATA c / .19730271018d2,     .209662681977d2,   -.483429455355d0,
     1         .605743189245d1,    .2256023885d2,     -.987532442d1,
     2        -.43135538513d1,     .458155781d0,      -.47754901883d-1,
     3         .41238460633d-2,   -.27929052852d-3,    .14481695261d-4,
     4        -.56473658748d-6,    .16200446d-7,      -.3303822796d-9,
     5         .451916067368d-11, -.370734122708d-13, 
     6         .137546068238d-15/


      tt  = t / 1.0d2
      tl  = DLOG(tt)
      gi  = -(c(1)/tt + c(2)) * tl
      hi  = (c(2) + c(1)*(1.0d0 - tl)/tt)
      cpi = c(2) - c(1)/tt

      DO 8 i=3,18
           emult = power(tt,DBLE(i-6))
           gi  = gi - c(i) * emult
           hi  = hi + c(i) * (i-6) * emult
 8         cpi = cpi + c(i) * (i-6) * (i-5) * emult

      ai  = gi - 1.0d0
      ui  = hi - 1.0d0
      cvi = cpi - 1.0d0
      si  = ui - ai

      RETURN
      END

******************************************************************************

*** dalHGK - Computes/returns (d(alpha)/dt)p(d,t,alpha) 
*            for the Haar et al. (1983) equation of state.
                
      DOUBLE PRECISION FUNCTION dalHGK(d,t,alpha)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      
      DOUBLE PRECISION tempi(4), densi(4), betai(4), alphai(4),
     1                 g(40), k, l, km, lm, kp, lp
      INTEGER          ll(40), kk(40)

      COMMON /aconst/ wm, gascon, tz, a, z, dz, y, uref, sref
      COMMON /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
      COMMON /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
      COMMON /resf/   ar, gr, sr, ur, hr, cvr, dpdtr
      COMMON /qqqq/   q, q5
      COMMON /nconst/ g, kk, ll, n
      COMMON /addcon/ tempi, densi, betai, alphai
      COMMON /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL

      SAVE


*** evaluate derivatives for the base function

      y     = .25d0 * b1 * d
      x     = 1.0d0 - y
      dydtp = (d/4.0d0)*(b1t - b1*alpha)

      dbdd = gascon*t * ((b1/4.0d0/x) * (1.0d0 - (g2-1.0d0)/x + 
     1       (g1+g2+1.0d0)/x/x) + b2 - b1*gf + 1.0d0/d)

      db2dd = gascon*t* ((b1*b1/16.0d0/x/x) * (1.0d0 - 
     1        2.0d0*(g2-1.0d0)/x + 3.0d0*(g1+g2+1.0d0)/x/x) - 
     2        1.0d0/d/d)

      db2ddt = gascon*t * ((b1t/4.0d0/x/x) *
     1         (1.0d0 - (g2-1.0d0)*(1.0d0+y)/x + 
     2         (g1+g2+1.0d0)*(1.0d0+2.0d0*y)/x/x) +
     3         b2t - gf*b1t) + dbdd/t 

      db2dtp = dbdd/t + gascon*t* ( (b1*dydtp/4.0d0/x/x/x) * 
     1         (1.0d0 - g2 + 2.0d0*(g1+g2+1.0d0)/x) + 
     2         ((x*b1t + b1*dydtp)/4.0d0/x/x) * 
     3         (1.0d0 - (g2-1.0d0)/x + (g1+g2+1.0d0)/x/x) + 
     4         b2t - gf*b1t + alpha/d )
 
      db3ddt = db2dd/t + gascon*t * ( (b1*b1*dydtp/8.0d0/x/x/x/x) * 
     1         (1.0d0 - g2 + 3.0d0*(g1+g2+1.0d0)/x) + 
     2         (b1*(x*b1t + b1*dydtp)/8.0d0/x/x/x) * 
     3         (1.0d0 - 2.0d0*(g2-1.0d0)/x + 3.0d0*(g1+g2+1.0d0)/x/x)
     4         - 2.0d0*alpha/d/d )
      
      db3dtt = (db2ddt - dbdd/t)/t + gascon*t* ( 
     1         (b1t*dydtp/2.0d0/x/x/x/x) * (1.0d0 - g2 + 
     2         (g1+g2+1.0d0)*(2.0d0+y)/x) + 
     3         ((x*b1tt + 2.0d0*b1t*dydtp)/4.0d0/x/x/x) * (1.0d0 - 
     4         (g2-1.0d0)*(1+y)/x + (g1+g2+1.0d0)*(1.0d0+2.0d0*y)/x/x)
     5         + b2tt - gf*b1tt ) + (t*db2dtp - dbdd)/t/t

***********************************************************

*** evaluate derivatives for the residual function

*      drdd   = q/d/d
*      dr2dd  = (q5 - 2.0d0/d*q)/d/d
*      dr2ddt = dpdtr/d/d

      e1  = DEXP(-a * d)
      e2  = 1.0d0 - e1
      tzt = tz / t

      drdd   = 0.0d0
      dr2dd  = 0.0d0
      dr2ddt = 0.0d0
      dr2dtp = 0.0d0
      dr3ddt = 0.0d0
      dr3dtt = 0.0d0

*** evaluate terms 1-36

      DO 10 i=1,n
           k = DBLE(kk(i)) + 1.0d0
           l = DBLE(ll(i)) - 1.0d0
           km = k - 1.0d0
           lm = l - 1.0d0
           kp = k + 1.0d0
           lp = l + 1.0d0
           xtzt = power(tzt,l)

           drdd   = drdd + g(i) * xtzt*power(e2,km)*e1

           dr2dd  = dr2dd + g(i) * e1*xtzt*power(e2,km) *
     1              (km*e1/e2 - 1.0d0)

           dr2ddt = dr2ddt - g(i)*e1*l*power(e2,km)*power(tzt,lp)/tz

           dr2dtp = dr2dtp + g(i)*e1*power(e2,km)*xtzt *
     1              ( d*alpha - l/t - km*e1*d*alpha/e2 )
          
           dr3ddt = dr3ddt + g(i)*( km*d*alpha*e1*e1*xtzt*
     1              power(e2,k-3.0d0) + e1*xtzt*power(e2,km)*
     2              (km*e1/e2 - 1.0d0) * (d*alpha -  l/t - 
     3              km*d*alpha*e1/e2) )

           dr3dtt = dr3dtt + g(i)*l*e1*power(e2,km)*power(tzt,lp)/tz
     1              * ( lp/t + d*alpha*km*e1/e2 - d*alpha )

 10        CONTINUE

*** evaluate terms 37-40

      DO 20 i=37,40
           k  = DBLE(kk(i))
           l  = DBLE(ll(i))
           km = k - 1.0d0
           lm = l - 1.0d0
           kp = k + 1.0d0
           lp = l + 1.0d0
           ai = alphai(i-36)
           bi = betai(i-36)
           di = densi(i-36)
           ti = tempi(i-36)
           tau = t/ti - 1.0d0
           del = d/di - 1.0d0
           IF (DABS(del) .LT. 1.0d-10)  del = 1.0d-10

           ex1 = -ai * power(del,k)
           IF (ex1 .LT. EXPTOL) THEN
                dex = 0.0d0
           ELSE
                dex = DEXP(ex1)
           END IF
           ex2  = -bi * tau * tau
           IF (ex2 .LE. EXPTOL) THEN
                tex = 0.0d0
           ELSE
                tex = DEXP(ex2)
           END IF
           ex12  = dex * tex
           qm    = l/del - k*ai*power(del,km)
           xdell = power(del,l)
           xdelk = power(del,k)

           drdd   = drdd + g(i)*xdell*ex12/di*qm

           dr2dd  = dr2dd + g(i)*xdell*ex12/di/di * (qm*qm - 
     1              l/di/di - ai*k*km*power(del,k-2.0d0))

           dr2ddt = dr2ddt - g(i)*2.0d0*bi*tau*ex12*xdell/ti/di*qm

           dr2dtp = dr2dtp + g(i)/di*( d*alpha*xdell*ex12/di/del/del *
     1              (l + ai*k*km*xdelk) + qm * ( ex12 *
     2              ( xdell* (k*ai*d*alpha*power(del,km)/di - 
     3              2.0d0*bi*tau/ti) - l*d*alpha*power(del,lm)/di) ) )
   
           dr3ddt = dr3ddt + g(i)/di/di*( xdell*ex12* (2.0d0*qm*
     1              (l*d*alpha/di/del/del + ai*k*km*d*alpha*
     2              power(del,k-2.0d0)/di) - 2.0d0*l*d*alpha/di/del
     3              /del/del + ai*k*km*(k-2.0d0)*power(del,k-3.0d0)*
     4              d*alpha/di) + (qm*qm - l/del/del - ai*k*km*
     5              power(del,k-2.0d0)) *(ex12*xdell*( ai*k*
     6              power(del,k-1.0d0)*d*alpha/di - 2.0d0*bi*tau/ti ) - 
     7              ex12*l*power(del,l-1.0d0)*d*alpha/di) )

           dr3dtt = dr3dtt - 2.0d0*g(i)*bi/ti/di * ( tau*xdell*ex12*d*
     1              alpha/del/del/di * (l + ai*k*km*power(del,k)) +
     2              qm*( xdell*ex12*( ai*k*d*alpha*tau*power(del,km)/di
     3              + (1.0d0 - 2.0d0*bi*tau*tau)/ti - 
     4              tau*l*d*alpha/di/del ) ) )


 20        CONTINUE

*** compute (d(alpha)/dT)P

      dalHGK = ((db3dtt + dr3dtt)*(2.0d0*(dbdd + drdd) + 
     1         d*(db2dd + dr2dd)) - 
     2         (db2ddt + dr2ddt)*(2.0d0*(db2dtp + dr2dtp) + 
     3         d*(db3ddt + dr3ddt) - d*alpha*(db2dd + dr2dd))) /
     4         (2.0d0*(dbdd + drdd) + d*(db2dd + dr2dd)) / 
     5         (2.0d0*(dbdd + drdd) + d*(db2dd + dr2dd))
     
      RETURN

      END

******************************************************************************

*** denHGK - Computes density (d in g/cm3) and dP/dD (dPdd) as
*            f(p(MPa),t(degK)) from an initial density guess (dguess).

      SUBROUTINE denHGK(d,p,dguess,t,dpdd)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      COMMON /qqqq/   q0, q5
      COMMON /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
      COMMON /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
      COMMON /RTcurr/ rt

      SAVE


      i  = 0
      d  = dguess

 10   i = i + 1

      IF (d .LE. 0.0d0)  d = 1.0d-8
      IF (d .GT. 1.9d0)  d = 1.9d0
  
      CALL resid(t,d)
      CALL base(d,t)

      pp   = rt * d * pb + q0
      dpdd = rt * (z + y * dz) + q5

*** if  dpdd < 0  assume d in 2-phase region and adjust accordingly ***

      IF (dpdd .GT. 0.0d0)        GO TO 20

      IF (dguess .GE. 0.2967d0)   d = d * 1.02d0
      IF (dguess .LT. 0.2967d0)   d = d * 0.98d0
      IF (i .LE. 10)              GO TO 10

 20   dpdx = dpdd * 1.1d0
      IF (dpdx .LT. 0.1d0)  dpdx = 0.1d0
      dp   = DABS(1.0d0 - pp/p)

      IF ((dp     .LT. 1.0d-8) .OR.
     1   ((dguess .GT. 0.3d0) .AND. (dp .LT. 1.0d-7)) .OR. 
     2   ((dguess .GT. 0.7d0) .AND. (dp .LT. 1.0d-6)))      RETURN

      x    = (p - pp) / dpdx
      IF (DABS(x) .GT. 0.1d0)  x = x * 0.1d0 / DABS(x)
      d = d + x
      IF (d .LE. 0.0d0)  d = 1.0d-8
      IF (i .LE. 30)    GO TO 10
      
      RETURN

      END

***********************************************************************

*** PsHGK - Returns an approximation to Psaturation(T) that agrees
*           to within 0.02% of that predicted by the HGK surface
*           for temperatures up to within roughly a degree of
*           the critical point.

      DOUBLE PRECISION FUNCTION PsHGK(t)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION  a(8)

      SAVE

      DATA a /-.78889166d1,  .25514255d1, -.6716169d1,  .33239495d2,
     1        -.10538479d3,  .17435319d3, -.14839348d3, .48631602d2/


      IF (T .LE. 314.0d0) THEN
           pl    = 6.3573118d0 - 8858.843d0/t + 
     1             607.56335d0 * power(t,-0.6d0)
           PsHGK = 0.1d0 * DEXP(pl)
      ELSE
           v = t / 647.25d0
           w = DABS(1.0d0 - v)
           b = 0.0d0
           DO 4 i=1,8
                z = i
                b = b + a(i)*power(w,(z + 1.0d0)/2.0d0)
 4              CONTINUE  
           q = b / v
           PsHGK = 22.093d0 * DEXP(q)
      END IF

      RETURN
      END

***********************************************************************

*** TsHGK - Returns Tsaturation(P).

      DOUBLE PRECISION FUNCTION TsHGK(p)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      SAVE


      TsHGK = 0.0d0

      IF (p .GT. 22.05d0)  RETURN

      k  = 0
      pl = 2.302585d0 + DLOG(p)
      tg = 372.83d0 + 
     1     pl*(27.7589d0 + pl*(2.3819d0 + pl*(0.24834d0 + 
     2     pl*0.0193855d0)))

 1    IF (tg .LT. 273.15d0)  tg = 273.15d0
      IF (tg .GT. 647.00d0)  tg = 647.00d0

      IF (k .GE. 8) THEN
           TsHGK = tg
      ELSE
           k  = k + 1
           pp = PsHGK(tg)
           dp = TdPsdT(tg)
           IF (ABS(1.0d0 - pp/p) .LT. 1.0d-5) THEN
                TsHGK = tg
           ELSE
                tg = tg * (1.0d0 + (p - pp)/dp)
                GO TO 1
           END IF
      END IF

      RETURN
      END

***********************************************************************

*** TdPsdT - Returns  T*(dPsat/dT).

      DOUBLE PRECISION FUNCTION TdPsdT(t)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION a(8)

      SAVE

      DATA a /-.78889166d1,  .25514255d1, -.6716169d1,  .33239495d2,
     1        -.10538479d3,  .17435319d3, -.14839348d3, .48631602d2/


      v = t / 647.25d0
      w = 1.0 - v
      b = 0.0d0
      c = 0.0d0

      DO 4 i=1,8
           z = i
           y = a(i) * power(w,(z + 1.0d0)/2.0d0)
           c = c + y/w*(0.5d0 - 0.5d0*z - 1.0d0/v)
 4         b = b + y

      q      = b / v
      TdPsdT = 22.093d0 * DEXP(q) * c

      RETURN

      END

***********************************************************************

*** corr - Computes liquid and vapor densities (dliq & dvap)
*          and  (Gl-Gv)/RT  (delg) for T-P conditions on or
*          near the saturation surface.

      SUBROUTINE corr(itripl,t,p,dl,dv,delg,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER epseqn

      COMMON /qqqq/   q00, q11
      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
      COMMON /fcts/   ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, dpdd,
     1                cjtt, cjth
      COMMON /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
      COMMON /RTcurr/ rt
      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
      COMMON /HGKcrt/ tcHGK, dcHGK, pcHGK

      SAVE


      CALL bb(t)

      dguess = dl
      IF (dl .LE. 0.0d0)  dguess = 1.11d0 - 0.0004d0*t

      CALL denHGK(dl,p,dguess,t,dpdd)
      CALL ideal(t)
      CALL thmHGK(dl,t)
*** save liquid properties
      CALL dimHGK(1,itripl,t,p,dl,epseqn)
      gl   = gd

      dguess = dv
      IF (dv .LE. 0.0d0)   dguess = p / rt

      CALL denHGK(dv,p,dguess,t,dpdd)
      IF (dv .LT. 5.0d-7)  dv = 5.0d-7
      CALL ideal(t)
      CALL thmHGK(dv,t)
*** vapor properties will be available
*** in COMMON /fcts/ (dimensionless) after
*** pcorr's final call of corr (delg < 10d-4)
      gv   = gd
      delg = gl - gv

      RETURN
      END

***********************************************************************

*** pcorr - Computes Psaturation(T) (p) and liquid and vapor
*           densities (dl & dv) from refinement of an initial
*           approximation (PsHGK(t)) in accord with  Gl = Gv.

      SUBROUTINE pcorr(itripl,t,p,dl,dv,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER epseqn

      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref

      SAVE

      p  = PsHGK(t)
      dl = 0.0d0
      dv = 0.0d0

 2    CALL corr(itripl,t,p,dl,dv,delg,epseqn)

      dp = delg * gascon * T / (1.0d0/dv - 1.0d0/dl)
      p  = p + dp
      IF (DABS(delg) .GT. 1.0d-4)  GO TO 2

      RETURN
      END

************************************************************

*** tcorr - Computes Tsaturation(P) (t) and liquid and vapor
*           densities (dl & dv) from refinement of an initial
*           approximation (TsHGK(p)) in accord with  Gl = Gv.

      SUBROUTINE tcorr(itripl,t,p,dl,dv,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER  epseqn

      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
      COMMON /RTcurr/ rt

      SAVE

      
      t = TsHGK(p)
      IF (t .EQ. 0.0d0) RETURN
      dl = 0.0d0
      dv = 0.0d0

 1    rt = t * gascon
      CALL corr(itripl,t,p,dl,dv,delg,epseqn)

      dp = delg * gascon * t / (1.0d0/dv - 1.0d0/dl)
      t = t * (1.0d0 - dp/TdPsdT(t))

      IF (DABS(delg) .GT. 1.0d-4) GO TO 1

      RETURN
      END

***************************************************************

*** LVSeqn - Computes thermodynamic and transport properties of 
*            critical region H2O (369.85-419.85 degC, 
*            0.20-0.42 gm/cm3) from the fundamental equation given 
*            by Levelt Sengers, et al (1983): J.Phys.Chem.Ref.Data,
*            V.12, No.1, pp.1-28.

      SUBROUTINE LVSeqn(isat,iopt,itripl,T,P,Dens,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP = 23)

      DOUBLE PRECISION  wprops(NPROP), wpliq(NPROP), Dens(2)
      LOGICAL           cpoint
      INTEGER           epseqn

      COMMON /coefs/  a(20), q(20), x(11)
      COMMON /crits/  Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /therm/  AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw,
     1                heat, Speed
      COMMON /satur/  Dliq, Dvap, DH2O, iphase
      COMMON /param/  r1, th1
      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
      COMMON /wpvals/ wprops, wpliq

      SAVE


      cpoint = .FALSE.
      DH2O = Dens(1)

 10   CALL LVSsat(iopt,isat,T,P,DH2O)

      IF ((isat .NE. 0) .OR. (iopt .NE. 1))  CALL denLVS(isat,T,P)

      IF (isat .EQ. 0) THEN
           Dens(1) = DH2O
      ELSE
           Dens(1) = Dliq
           Dens(2) = Dvap
      END IF
                 
      IF (isat .EQ. 0) THEN
           CALL thmLVS(isat,T,r1,th1)
           CALL dimLVS(isat,itripl,th1,T,P*1.0d1,dl,dv,wprops,epseqn)
           IF (cpoint) THEN
                CALL cpswap
                Dens(1) = cdens
                Dens(2) = cdens
                isat = 1
                iopt = ioptsv
           END IF
      ELSE
           th1 = -1.0d0
           CALL thmLVS(isat,T,r1,th1)
           CALL dimLVS(isat,itripl,th1,T,P*1.0d1,dl,dv,wprops,epseqn)
           th1 =  1.0d0
           CALL thmLVS(isat,T,r1,th1)
           CALL dimLVS(isat,itripl,th1,T,P*1.0d1,dl,dv,wpliq,epseqn)
           IF (dl .EQ. dv) THEN
                cpoint = .TRUE.
                cdens = dl
                T = 647.0670000003d0
                P =  22.0460000008d0
                ioptsv = iopt
                iopt = 2
                isat = 0
                GO TO 10 
           END IF
      END IF

      END

*********************************************************************

*** cpswap - Load critical point A, G, U, H, S, Vs, Di, ZB, 
*            albe values from wpliq into wprops and 
*            approximations to critical Cv, Cp, alpha, beta, 
*            visc, tcond, Prndtl, tdiff, visck, YB, QB, XB, 
*            daldT, st values from wprops into wpliq.

      SUBROUTINE cpswap

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      
      PARAMETER (NPROP = 23)

      INTEGER           aw, gw, sw, uw, hw, cvw, cpw, vsw, alw, bew,
     1                  diw, viw, tcw, stw, tdw, Prw, vikw, albew,
     2                  ZBw, YBw, QBw, dalwdT, XBw
      DOUBLE PRECISION  wprops(NPROP), wpliq(NPROP)

      COMMON /wpvals/ wprops, wpliq
      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc

      SAVE

      DATA aw, gw, sw, uw, hw, cvw, cpw, vsw, alw, bew, diw, viw,
     1     tcw, stw, tdw, Prw, vikw, albew, ZBw, YBw, QBw, 
     2     dalwdT, XBw
     2   /  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
     3     13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 /

      
      wprops(aw)    = wpliq(aw)
      wprops(gw)    = wpliq(gw)
      wprops(sw)    = wpliq(sw)
      wprops(uw)    = wpliq(uw)
      wprops(hw)    = wpliq(hw)
      wprops(diw)   = wpliq(diw)
      wprops(ZBw)   = wpliq(ZBw)
      wprops(stw)   = wpliq(stw)

      wpliq(cvw)    = wprops(cvw)
      wpliq(cpw)    = wprops(cpw)
      wpliq(alw)    = wprops(alw)
      wpliq(bew)    = wprops(bew)
      wpliq(YBw)    = wprops(YBw)
      wpliq(QBw)    = wprops(QBw)
      wpliq(XBw)    = wprops(XBw)
      wpliq(tcw)    = wprops(tcw)
      wpliq(tdw)    = wprops(tdw)
      wpliq(Prw)    = wprops(Prw)
      wpliq(dalwdT) = wprops(dalwdT)
      wpliq(albew)  = wprops(albew)

      wprops(vsw)   = 0.429352766443498d2 * fs
      wprops(viw)   = 1.0d6
      wprops(vikw)  = 1.0d6

      wpliq(vsw)    = wprops(vsw)
      wpliq(viw)    = wprops(viw)
      wpliq(vikw)   = wprops(vikw)

      END      

*********************************************************************

*** LVSsat - If  isat=1,  computes  Psat(T) or Tsat(P) (iopt=1,2).
*            If  isat=0,  checks whether  T-D or T-P (iopt=1,2)
*            falls on or within  TOL  of the liq-vap surface; if so,
*            isat <- 1  and  T <- Tsat.

      SUBROUTINE LVSsat(iopt,isat,T,P,D)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      COMMON /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
      COMMON /crits/  Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon

      SAVE

      DATA ERRTOL, TCTOL / 1.0d-12, 1.0d-2 /


      IF (isat .EQ. 1) THEN
           IF (iopt .EQ. 1) THEN
                P = Pfind(isat,T,D)
           END IF
           T = TsLVS(isat,P)
      ELSE
           IF (iopt .EQ. 1) THEN
                P = Pfind(isat,T,D)
           END IF
           IF (P-ERRTOL .GT. Pc) THEN
                RETURN
           ELSE
                CALL backup
                Tsat = TsLVS(isat,P)
                IF (DABS(Tsat-T) .LT. TCTOL) THEN
                     T = Tsat
                     isat = 1
                ELSE
                     CALL restor
                END IF
           END IF
      END IF

      RETURN
      END

*********************************************************************

*** denLVS - Calculates  DH2O(T,P)  or  Dvap,Dliq(T,P) from the 
*            Levelt Sengers, et al (1983) critical region 
*            equation of state.

      SUBROUTINE denLVS(isat,T,P)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION s(2), sd(2)

      COMMON /coefs/ a(20), q(20), x(11)
      COMMON /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /satur/ Dliq, Dvap, DH2O, iphase
      COMMON /therm/ AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw,
     1               heat, Speed
      COMMON /param/ r1, th1
      COMMON /deri2/ dPdD, dPdT

      SAVE

      EQUIVALENCE (Dmin, x(4)), (Dmax, x(5)), (pw11, q(9)),  
     1            (xk0,  a(7)), (xk1,  a(12))
      
      IF (isat .EQ. 0) THEN
           DH2O = rhoc
           DO 10 i=1,20
                Pnext = Pfind(isat,T,DH2O)
                Pdif  = Pnext - P
                IF (iphase .EQ. 2) THEN
                     IF (DABS(Pdif) .LE. 0.0d0) THEN
                          RETURN
                     ELSE
                     END IF
                     IF (Pdif .LT. 0.0d0) THEN
                          DH2O = Dmax
                     ELSE
                          DH2O = Dmin
                     END IF
                ELSE
                     delD  = -Pdif/dPdD
                     DH2O = DH2O + delD
                     IF (DH2O .LT. Dmin)  DH2O = Dmin
                     IF (DH2O .GT. Dmax)  DH2O = Dmax
                     IF (DABS(delD/DH2O) .LT. 1.0d-6)  RETURN
                END IF
 10        CONTINUE
      ELSE
           Tw   = -Tc/T
           dTw  = 1.0d0 + Tw

           CALL ss(r1,th1,s,sd)
           rho1 = 1.0d0+pw11*dTw+a(1)*(s(1)+s(2))
           rho2 = xk0*power(r1,a(6)) + xk1*power(r1,q(16))

           Dvap = rhoc * (rho1 - rho2)
           Dliq = rhoc * (rho1 + rho2)

           RETURN
      END IF

      RETURN
      END

*********************************************************************

*** TsLVS - Returns saturation T(P) 

      DOUBLE PRECISION FUNCTION TsLVS(isat,P)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      COMMON /therm/ AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw,
     1               heat, Speed
      COMMON /satur/ Dliq, Dvap, DH2O, iphase
      COMMON /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /deri2/ dPdD, dPdT

      SAVE

      
      TsLVS2 = Tc - 1.0d0
      D = rhoc

      DO 10 i=1,20
           Pnext = Pfind(isat,TsLVS2,D)
           dT = (Pnext - P)/dPdT
           TsLVS2 = TsLVS2 - dT
           IF (TsLVS2 .GT. Tc) THEN
                TsLVS2 = Tc
           ELSE
                IF (DABS(dT/TsLVS2) .LT. 1.0d-8) THEN
                     GO TO 20
                ELSE
                END IF
           END IF
 10   CONTINUE

 20   TsLVS = TsLVS2

      RETURN
      END

*********************************************************************

*** Pfind - Returns P(T,D).  Computes (dP/dD)T when invoked by SUB 
*           Dfind (isat=0) and (dP/dT)D when invoked by SUB TsLVS 
*           (isat=1).  Also computes 1st & 2nd partial derivatives 
*           the singular part of the potential (Delta P tilde) that
*           are used in SUB thmLVS.

      DOUBLE PRECISION FUNCTION Pfind(isat,T,D)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION s(2), xk(2), sd(2)

      COMMON /coefs/ a(20), q(20), x(11)
      COMMON /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /satur/ Dliq, Dvap, DH2O, iphase
      COMMON /therm/ AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw,
     1               heat, Speed
      COMMON /param/ r1, th1
      COMMON /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
      COMMON /deriv/ amu, s, Pw, Tw, dTw, dM0dT, dP0dT,
     1               d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
      COMMON /deri2/ dPdD, dPdT
***************************************
      COMMON /abc2/  r, th
***************************************

      SAVE

      EQUIVALENCE (Pw1, a(5)),   (Pw2, a(4)),   (Pw3,  a(2)),
     1            (amc, a(13)),  (am1, a(14)),  (am2,  a(15)),
     2            (am3, a(16)),  (p00, q(11)),  (p20,  q(12)),
     3            (p40, q(13)),  (p01, q(18)),  (p21,  q(19)),
     4            (p41, q(20)),  (aa,  a(10)),  (xk0,  a(7)),
     5            (xk1, a(12)),  (pw11,q(9)),   (alpha,q(10)),
     6            (alhi,q(15)),  (besq,a(9))

      xk(1) = xk0
      xk(2) = xk1
      IF (DABS(T-Tc) .LT. FPTOL)  T = Tc
      Tee   = (T-Tc)/Tc
      Tw    = -Tc/T
      dTw   = 1.0d0 + Tw

      IF (isat .EQ. 0) THEN
           rho = D / rhoc
           CALL conver(rho,Tee,amu,th1,r1,rho1,s,rhodi,err)
      ELSE
           th1 = -1.0d0
           th  = th1
           r1  = dTw/(1.0d0-besq)
           r   = r1
           CALL ss(r1,th1,s,sd)
           rho = th1 * (xk0*power(r1,a(6)) + 
     1           xk1*power(r1,q(16))) +
     2           a(1)*(s(1)+s(2))
           rho = 1.0d0+pw11*dTw+rho
           amu = 0.0d0
           D = rho * rhoc
      END IF

      tt1 = th1*th1
      tt2 = tt1*tt1

      Pw0  = 1.0d0+dTw*(Pw1+dTw*(Pw2+dTw*Pw3))

      IF (isat .EQ. 0) THEN
           Pwmu = amu*rhodi
      ELSE
           Pwmu = 0.0d0
      END IF

      p0th = p00+p20*tt1+p40*tt2
      p1th = p01+p21*tt1+p41*tt2

      dPw0 = xk0*p0th*power(r1,2.0d0-alpha)
      dPw1 = xk1*p1th*power(r1,2.0d0-alhi)
      dPw  = aa*(dPw0+dPw1)

      Pw   = Pw0 + Pwmu + dPw

      Pfind = Pw * Pcon * T

      IF (DABS(th1) .LT. 1.0d0) THEN
           iphase = 1
      ELSE
           iphase = 2

           dP0dT = Pw1+dTw*(2.0d0*Pw2+3.0d0*Pw3*dTw)
           dM0dT = am1+dTw*(2.0d0*am2+3.0d0*am3*dTw)
           Uw    = dP0dT-rho*dM0dT+pw11*amu+s(1)+s(2)

           dPdTcd = Uw + rho*dM0dT
           dPwdTw = Pw - Tw*dPdTcd

           dPdT   = Pcon * dPwdTw

      END IF

      CALL aux(r1,th1,d2PdT2,d2PdMT,d2PdM2,aa,xk,sd,Cvcoex)

      IF (iphase .EQ. 1) dPdD = dPcon * D * T / d2PdM2

      RETURN
      END

***************************************************************

*** aux - Calculates some second derivatives of the 
*         anomalous part of the equation of state.

      SUBROUTINE aux(r1,th1,d2PdT2,d2PdMT,d2PdM2,aa,xk,sd,Cvcoex)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION xk(2), s(2), sd(2), w(2), y(2), z(2), coex(2)

      COMMON /coefs/ a(20), q(20), x(11)

      SAVE

      EQUIVALENCE (cc,   a(1)),  (beta, a(6)),  (besq,a(9)),
     1            (delta,a(11)), (alpha,q(10)), (s00, a(17)),
     2            (s20,  a(18)), (s01,  a(19)), (s21, a(20))


      deli  = 0.0d0
      s(1)   = s00+s20*th1*th1
      s(2)   = s01+s21*th1*th1
      sd(1)  = 2.0*th1*s20
      sd(2)  = 2.0*th1*s21
      ww     = 0.0d0
      yy     = 0.0d0
      zz     = 0.0d0
      gamma  = beta*(delta-1.0d0)
      tt1    = th1*th1
      ter    = 2.0d0*beta*delta-1.0d0
      g      = (1.0+(besq*ter-3.0)*tt1 - besq*(ter-2.0)*tt1*tt1)
      Cvcoex = 0.0d0

      DO 30 i=1,2
           alhi    = alpha - deli
           beti    = beta + deli
           gami    = gamma - deli
           IF (r1 .NE. 0.0d0) THEN
                w(i)    = (1.0-alhi)*(1.0-3.0*tt1)*s(i) - 
     1                    beta*delta*(1.0-tt1)*th1*sd(i)
                w(i)    = (w(i)*power(r1,-alhi))/g
                w(i)    = w(i) * xk(i)
                ww      = ww + w(i)

                y(i)    = beti*(1.0d0-3.0d0*tt1)*th1 - 
     1                    beta*delta*(1.0d0-tt1)*th1
                y(i)    = (y(i)*power(r1,beti-1.0d0)) * xk(i) / g
                yy      = yy + y(i)

                z(i)    = 1.0d0-besq*(1.0d0-(2.0d0*beti))*tt1
                z(i)    = (z(i)*power(r1,-gami)) * xk(i) / g
                zz      = zz + z(i)

                a1 = (beta*(delta-3.0d0)-3.0d0*deli-besq*alhi*gami) / 
     1               (2.0d0*besq*besq*(2.0d0-alhi)*(1.0d0-alhi)*alhi)
                a2 = 1+((beta*(delta-3.0d0)-3.0d0*deli-besq*alhi*ter) /
     1                  (2.0d0*besq*(1.0d0-alhi)*alhi))
                a2 = -a2
       
                a4 = 1.0d0+((ter-2.0d0)/(2.0d0*alhi))
                f1 = a1 + a2 + a4
           
                coex(i) = ((2.0d0-alhi)*(1.0d0-alhi)*power(r1,-alhi) * 
     1                    f1*xk(i))
                Cvcoex  = Cvcoex + coex(i)
           END IF
           deli    = 0.5d0
 30   CONTINUE

      d2PdT2 = aa * ww
      d2PdMT = yy + aa*cc*ww
      d2PdM2 = zz/aa + 2.0d0*cc*yy + cc*cc*aa*ww

      RETURN
      END

***************************************************************

*** conver - Transforms  T,D  to  parametric variables  r,theta 
*            according to the revised and scaled equations.

      SUBROUTINE conver(rho,Tee,amu,th1,r1,rho1s,s1,rhodi,error1)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
    
      DOUBLE PRECISION s1(2), sd(2)

      COMMON /coefs/ a(20), q(20), x(11)
      COMMON /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
**************************************************************
      COMMON /abc2/  r, th
**************************************************************

      SAVE

      EQUIVALENCE (beta,a(6)),  (delta,a(11)),  (xk1,  a(12)),
     1            (cc,  a(1)),  (alhi, q(15)),  (alpha,q(10)),
     2            (besq,a(9)),  (p11,  q(9)),   (deli, q(14)),
     3            (p1w, q(18)), (p2w,  q(19)),  (p4w,  q(20)),
     4            (aa,  a(10)), (xk0,  a(7)),   (s00,  a(17)),
     5            (s20, a(18)), (betai,q(16))


      Tstar  = Tee + 1.0d0
      dtstin = 1.0d0 - (1.0d0 / Tstar)
      r1     = dtstin
     
      IF (dtstin .LE. 0.0d0)  THEN
           r1  = dtstin/(1.0d0-besq)
           th1 = 1.0d0
      ELSE
           th1 = 0.0d0
      END IF

      CALL ss(r1,th1,s1,sd)

      rhodi  = 1.0d0 + p11*dtstin
      rhodit = rhodi + cc*s1(1) + cc*s1(2)
      drho   = rho - rhodit
      amu    = 0.0d0

      IF (dtstin .LE. 0.0d0) THEN
           rho1co = xk0*power(r1,beta) + xk1*power(r1,betai)
           twofaz = rho1co
           IF (DABS(drho) .LE. twofaz) THEN
                rho1s  = DSIGN(rho1co,drho) + cc*s1(1)
                th1    = DSIGN(1.00d0,drho)
                error1 = 1.0d0
                r = r1
                th = th1
                RETURN
           END IF
      END IF

      IF (drho .EQ. 0.0d0) THEN
           th1   = 0.0d0
           r1    = dtstin
           rho1s = cc*s1(1)
      END IF 

*** rule for first pass ***

      y1   = dtstin      
      den1 = rho - rhodit
 
      CALL rtheta(r1,th1,den1,y1)

      tt   = th1*th1
      amu  = aa*power(r1,beta*delta)*th1*(1.0d0-tt)
      y1   = dtstin + cc*amu

      CALL ss(r1,th1,s1,sd)

      rhoweg = xk1*power(r1,betai)*th1 + cc*s1(2)      
      rho1s  = den1 + cc*s1(1) + rhoweg
      error1 = rho - rhodi - rho1s
      r  = r1
      th = th1

      IF (DABS(error1) .LT. 1.0d-5) THEN
           RETURN
      END IF

*** rule for second pass ***

      den12 = rho - rhodi - cc*s1(1) + rhoweg
      
      IF (den12 .EQ. den1) den12 = den1 - 1.0d-6
   
      CALL rtheta(r1,th1,den12,y1)

      tt  = th1*th1
      amu = aa*power(r1,beta*delta)*th1*(1.0d0-tt)
      y1  = dtstin + cc*amu

      CALL ss(r1,th1,s1,sd)

      rhoweg = xk1*power(r1,betai)*th1 + cc*s1(2)
      rho1s2 = den12 + cc*s1(1) + rhoweg
      error2 = rho - rhodi - rho1s2

      IF (DABS(error2) .LE. 1.0d-5) THEN
           r  = r1
           th = th1
           error1 = error2
           rho1s  = rho1s2
           RETURN      
      END IF

*** rule for nth pass ***

      den2   = den12
      
      DO 44 isig=1,10
           slope  = (error2-error1)/(den2-den1)
           hold   = den2
           den2   = den1 - (error1/slope)
           
           CALL rtheta(r1,th1,den2,y1)

           tt  = th1*th1
           amu = aa*power(r1,beta*delta)*th1*(1.0d0-tt)
           y1  = dtstin + cc*amu

           CALL ss(r1,th1,s1,sd)

           rhoweg = xk1*power(r1,betai)*th1 + cc*s1(2)
           rho1s  = den2 + cc*s1(1) + rhoweg
           error1 = error2
           error2 = rho - rhodi - rho1s
           r  = r1
           th = th1

           IF (DABS(error2) .LT. 1.0d-6) RETURN
  
           den1 = hold

 44   CONTINUE

      RETURN
      END

*********************************************************************

*** rtheta - Fits data for  1.0 < theta < 1.000001.
*            Solves:
*                     rho = em*theta*(r**beta)
*                     Tee = r*(1.0d0-besq*theta*theta)
*
*   Routine given by Moldover (1978): Jour. Res. NBS, v. 84, n. 4,
*   p. 329 - 334.


      SUBROUTINE rtheta(r,theta,rho,Tee)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      COMMON /coefs/ a(20), q(20), x(11)

      SAVE

      EQUIVALENCE (beta,a(6)), (em,a(7)), (besq,a(9))


      IF (em .LE. 0.0d0  .OR.  besq .LE. 1.0d0) GO TO 600

      absrho = DABS(rho)

      IF (absrho .LT. 1.0d-12) GO TO 600

      bee = DSQRT(besq)

      IF (DABS(Tee) .LT. 1.0d-12) GO TO 495
      IF (Tee .LT. 0.0d0) THEN
           z = 1.0d0-(1.0d0-bee)*Tee/(1.0d0-besq) * 
     1         power(em/absrho,1.0d0/beta)
      ELSE
           z = power(1.0d0+Tee*power(em/bee/absrho,1.0d0/beta),
     1         -beta)
      END IF
      IF (z .GT. 1.00234d0*bee) GO TO 496

      c = -rho*bee/em/power(DABS(Tee),beta)
      z = DSIGN(z,rho)

      DO 500 n=1,16
           z2 = z*z
           z3 = 1.0d0 - z2
           dz = z3*(z+c*power(DABS(z3),beta))/(z3+2.0d0*beta*z2)
           z  = z - dz

           IF (DABS(dz/z) .LT. 1.0d-12) GO TO 498

 500  CONTINUE

 601  IF (DABS(theta) .GT. 1.0001d0) theta = theta/DABS(theta)
      RETURN

 498  theta = z/bee
      r     = Tee/(1.0d0-z*z)
      r     = DABS(r)
      RETURN

 495  theta = DSIGN(1.0d0,rho)/bee
      r     = power(rho/(em*theta),1.0d0/beta)
      RETURN

 496  theta = DSIGN(1.0d0,rho)
      r     = Tee/(1.0d0-besq)
      r     = DABS(r)
      RETURN

 600  IF (DABS(Tee) .LT. 1.0d-12) GO TO 601

      IF (Tee .LT. 0.0d0) GO TO 496

      theta = 1.0d-12
      r     = Tee
      RETURN

      END

*********************************************************************

*** ss - Computes terms of the summation that defines  dPotl/dT
*        and the 1st derivative of the theta (s) square polynomial.

      SUBROUTINE ss(r,th,s,sd)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION s(2), sd(2), sx(2)
    
      COMMON /coefs/ a(20), q(20), x(11)
***************************************************************
      COMMON /abc1/  dPdM
***************************************************************

      SAVE

      EQUIVALENCE (alpha,q(10)),  (beta,a(6)),  (besq,a(9)),
     1            (delta,a(11)),  (deli,q(14)), (alhi,q(15)),
     2            (beti, q(16)),  (gami,q(17)), (p00, q(11)),
     3            (p01,  q(18)),  (s00, a(17)), (s20, a(18)),
     4            (s01,  a(19)),  (s21,  a(20))

      tt    = th*th
      sx(1)  = s00 + s20*tt
      sd(1) = 2.0d0*s20*th
      sx(2)  = s01 + s21*tt
      sd(2) = 2.0d0*s21*th
      s(1)  = sx(1)*a(10)*a(7)*power(r,1.0d0-alpha)
      s(2)  = sx(2)*a(10)*a(12)*power(r,1.0d0-alhi)

      dPdM  = power(r,beta)*a(7)*th  + a(1)*power(r,1.0d0-alpha)*
     1        a(10)*a(7)*sx(1) +
     2        power(r,beti)*a(12)*th + a(1)*power(r,1.0d0-alhi)*
     3        a(10)*a(12)*sx(2)

      RETURN
      END

*****************************************************************

*** thmLVS - Calculates thermodynamic and transport properties
*            of critical region H2O using the Levelt Sengers, et al
*            (1983) equation of state.

      SUBROUTINE thmLVS(isat,T,r1,th1)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION s(2), xk(2), sd(2)

      COMMON /coefs/ a(20), q(20), x(11)
      COMMON /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /therm/ AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw,
     1               heat, Speed
      COMMON /satur/ Dliq, Dvap, DH2O, iphase
      COMMON /deriv/ amu, s, Pw, Tw, dTw, dM0dT, dP0dT,
     1               d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
      COMMON /deri2/ dPdD, dPdT
*************************************************************
      COMMON /abc1/  dPdM
      COMMON /abc3/  dPdTcd
*************************************************************

      SAVE

      EQUIVALENCE (pw2, a(4)),   (pw3, a(2)),  (besq,  a(9)),
     1            (amc, a(13)),  (am1, a(14)), (am2,   a(15)),
     2            (aa,  a(10)),  (xk0, a(7)),  (am3,   a(16)),
     3            (xk1, a(12)),  (pw11,q(9)),  (alpha, q(10)),
     4            (alhi,q(15)),  (pw1, a(5))

      d2P0dT = 2.0d0*pw2 + 6.0d0*pw3*dTw
      d2M0dT = 2.0d0*am2 + 6.0d0*am3*dTw

      dP0dT  = pw1+dTw*(2.0d0*pw2+3.0d0*pw3*dTw)
      dM0dT  = am1+dTw*(2.0d0*am2+3.0d0*am3*dTw)

      IF (isat .EQ. 0) THEN
           rho    = DH2O / rhoc
           Uw     = dP0dT-rho*dM0dT+pw11*amu+s(1)+s(2)
      ELSE
           rho    = th1 * (xk0*power(r1,a(6)) + xk1*power(r1,q(16)))
     1              + a(1)*(s(1)+s(2))
           rho    = 1.0d0+pw11*dTw+rho
           Uw     = dP0dT-rho*dM0dT+pw11*amu+s(1)+s(2)
           DH2O   = rho * rhoc
           dPdT2  = Pw - Tw*(Uw+rho*dM0dT)
           heat   = 1.0d3*T*(Pcon*dPdT2)*(1.0d0/Dvap-1.0d0/Dliq)

           CALL ss(r1,th1,s,sd)
           CALL aux(r1,th1,d2PdT2,d2PdMT,d2PdM2,aa,xk,sd,Cvcoex)
           IF (r1 .NE. 0.0d0) THEN
                dPdD = dPcon * DH2O * T / d2PdM2
           END IF
      END IF

      IF (r1 .NE. 0.0d0) THEN
           dPdTcd = dP0dT+pw11*(amu-rho/d2PdM2)+s(1)+s(2) - 
     1              d2PdMT*rho/d2PdM2
           dPwdTw = Pw - Tw*dPdTcd
           dPdTal = Pcon * dPwdTw

           CviTw2 = d2P0dT - rho*d2M0dT + d2PdT2 - 
     1              (pw11+d2PdMT)*(pw11+d2PdMT)/d2PdM2
           Cvw    = CviTw2 * Tw*Tw
           Cpw    = Cvw + d2PdM2*dPwdTw*dPwdTw / (rho*rho)
           betaw  = 1.0d0 / (DH2O*dPdD)
           alphw  = betaw * dPdTal
           Speed   = 1.0d3 * DSQRT(Cpw/Cvw*dPdD)
      ELSE
           Cvw   = 1.0d0
           Cpw   = 1.0d0
           betaw = 1.0d0
           alphw = 1.0d0
           Speed = 0.0d0
      END IF

      Hw = Pw - Tw*Uw
      Sw = Hw - rho*(amu+amc+dTw*(am1+dTw*(am2+dTw*am3)))

      Scond  = Scon/DH2O

      U      = Uw * Ucon/DH2O
      H      = Hw * Scond * T
      entrop = Sw * Scond
      AE     = U - T * entrop
      GE     = H - T * entrop
      Cv     = Cvw * Scond
      Cp     = Cpw * Scond

      RETURN
      END

********************************************************

*** dalLVS - Computes/returns (d(alpha)/dt)p(D,T,alpha) 
*            for the Levelt Sengers et al. (1983) 
*            equation of state.  Note that D (kg/m**3),
*            T (degK), P (MPa), alpha (degK**-1).

                
      DOUBLE PRECISION FUNCTION dalLVS(D,T,P,alpha)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION sss(2), xk(2), s(2), dsdT(2), sp(2), dspdT(2),
     1                 k(2), calpha(2), cbeta(2), cgamma(2),
     2                 u(2), v(2), w(2), dudT(2), dvdT(2), dwdT(2)

      COMMON /coefs/ aa(20), qq(20), xx(11)
      COMMON /crits/ Tc, Dc, Pc, Pcon, Ucon, Scon, dPcon
      COMMON /deriv/ amu, sss, Pw, Tw, dTw, dM0dT, dP0dT,
     1               d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
      COMMON /deri2/ dPdD, dPdT
*************************************************************
      COMMON /abc1/  dPdM
      COMMON /abc2/  r,th
      COMMON /abc3/  dPdTcd
*************************************************************

      SAVE

      EQUIVALENCE (a,   aa(10)), (c,   aa(1)),  (delta,  aa(11)), 
     1            (bsq, aa(9)),  (P11, qq(9)),  (Delta1, qq(14)),
     2            (P1,  aa(5)),  (P2,  aa(4)),  (P3,     aa(2)),
     3            (s00, aa(17)), (s01, aa(19)), (s20,    aa(18)), 
     4            (s21, aa(20))
 

      IF (r .EQ. 0.0d0) THEN
           dalLVS = 1.0d6
           RETURN
      END IF

      k(1)      = aa(7)
      k(2)      = aa(12)
      calpha(1) = qq(10)
      calpha(2) = qq(15)
      cbeta(1)  = aa(6)
      cbeta(2)  = qq(16)
      cgamma(1) = cbeta(1)*(delta - 1.0d0)
      cgamma(2) = cgamma(1) - Delta1
      delT      = (T - Tc) / T

      s(1)      = s00 + s20*th**2
      s(2)      = s01 + s21*th**2
      sp(1)     = 2.0d0*s20*th
      sp(2)     = 2.0d0*s21*th
      
*********************************************************************
***
*** Compute drdT and d0dT from solution of the linear system
***
***                      ax = b
***
*** d(dPdM)/dT = -D/Dc*alpha - P11*Tc/T**2 = ar1*drdT + a01*d0dT = b1
*** d(delT)/dT =           Tc/T**2         = ar2*drdT + a02*d0dT = b2
***

      b1 = -D/Dc*alpha - P11*Tc/T/T
      b2 =  Tc/T**2   

      ar1 = 0.0d0                 
      a01 = 0.0d0
      DO 10 i = 1,2
           ar1 = ar1 + k(i) * (cbeta(i)*th*power(r,cbeta(i)-1.0d0) + 
     1           a*c*(1.0d0 - calpha(i))*power(r,-calpha(i))*s(i))
           a01 = a01 + k(i) * (power(r,cbeta(i)) + a*c*sp(i)*
     1           power(r,1.0d0-calpha(i)))
 10        CONTINUE
      
      ar2 = 1.0d0 - bsq*th**2 - a*c*cbeta(1)*delta*
     1      (1.0d0 - th**2)*th*power(r,(cbeta(1)*delta - 1.0d0))
      a02 = 3.0d0*a*c*th**2*power(r,cbeta(1)*delta) - 
     1      2.0d0*bsq*r*th - a*c*power(r,cbeta(1)*delta)

*********************************************************************
*** solve the linear system with simplistic GE w/ partial pivoting
*********************************************************************

      IF (DABS(ar1) .GT. DABS(ar2)) THEN 
           amult = -ar2 / ar1
           d0dT  = (b2 + amult*b1) / (a02 + amult*a01)
           drdT  = (b1 - a01*d0dT) / ar1
      ELSE
           amult = -ar1 / ar2
           d0dT  = (b1 + amult*b2) / (a01 + amult*a02)
           drdT  = (b2 - a02*d0dT) / ar2
      END IF

*********************************************************************
***
*** Compute theta polynomials and their tempertaure derivatives
***

      dsdT(1)   = 2.0d0*s20*th*d0dT
      dsdT(2)   = 2.0d0*s21*th*d0dT
      dspdT(1)  = 2.0d0*s20*d0dT
      dspdT(2)  = 2.0d0*s21*d0dT

      q     = 1.0d0 + (bsq*(2.0d0*cbeta(1)*delta - 1.0d0) - 3.0d0)*
     1        th**2 - bsq*(2.0d0*cbeta(1)*delta - 3.0d0)*th**4

      dqdT  = 2.0d0*(bsq*(2.0d0*cbeta(1)*delta - 1.0d0) - 3.0d0)*
     1        th*d0dT - 4.0d0*bsq*(2.0d0*cbeta(1)*delta - 3.0d0)*
     2        th**3*d0dT

      DO 20 i = 1,2
           u(i)    = (1.0d0 - bsq*(1.0d0 - 2.0d0*cbeta(i))*th**2) / q
           dudT(i) = (-2.0d0*bsq*(1.0d0 - 2.0d0*cbeta(i))*th*d0dT - 
     1               u(i)*dqdT) / q
           v(i)    = ((cbeta(i) - cbeta(1)*delta)*th + 
     1               (cbeta(1)*delta - 3.0d0*cbeta(i))*th**3) / q
           dvdT(i) = ((cbeta(i) - cbeta(1)*delta)*d0dT + 
     1               3.0d0*(cbeta(1)*delta - 3.0d0*cbeta(i))*
     2               th**2*d0dT - v(i)*dqdT) / q
           w(i)    = ((1.0d0 - calpha(i))*(1.0d0 - 3.0d0*th**2)*
     1               s(i) - cbeta(1)*delta*(th - th**3)*sp(i)) / q
           dwdT(i) = ((1.0d0 - calpha(i))*((1.0d0 - 3.0d0*th**2)*
     1               dsdT(i) - 6.0d0*th*s(i)*d0dT) - cbeta(1)*
     2               delta*((th - th**3)*dspdT(i) + sp(i)*
     3               (d0dT - 3.0d0*th**2*d0dT)) - w(i)*dqdT) / q
 20        CONTINUE

*********************************************************************
***
*** Compute dP0dTT, ddelMT, dPdTT, dPdMMT, dPdMTT, dPPTT
***

      dP0dTT = Tc/T**2 * (2.0d0*P2 + 6.0d0*P3*delT) 

      ddelMT = a*power(r,cbeta(1)*delta)* (cbeta(1)*delta*th/r*
     1         (1.0d0 - th**2)*drdT + (1.0d0 - 3.0d0*th**2)*d0dT)

      dPdTT  = 0.0d0
      dPdMMT = 0.0d0
      dPdMTT = 0.0d0
      DO 30 i = 1,2
           dPdTT  = dPdTT + a*k(i) * (power(r,1.0d0-calpha(i))*
     1              dsdT(i) + s(i)*(1.0d0 - calpha(i))*
     2              power(r,-calpha(i))*drdT)
        
           dPdMMT = dPdMMT + k(i) * ((power(r,-cgamma(i))*dudT(i) -
     1              u(i)*cgamma(i)*power(r,-1.0d0-cgamma(i))*drdT) /
     2              a + 2.0d0*c*(power(r,cbeta(i)-1.0d0)*dvdT(i) + 
     3              v(i)*(cbeta(i) - 1.0d0)*power(r,cbeta(i)-2.0d0)*
     4              drdT) + a*c**2*(power(r,-calpha(i))*dwdT(i) - 
     5              calpha(i)*w(i)*power(r,-1.0d0-calpha(i))*drdT))

           dPdMTT = dPdMTT + k(i) * (power(r,cbeta(i)-1.0d0)*dvdT(i) +
     1              v(i)*(cbeta(i) - 1.0d0)*power(r,cbeta(i)-2.0d0)*
     2              drdT + a*c*(power(r,-calpha(i))*dwdT(i) - 
     3              calpha(i)*power(r,-1.0d0-calpha(i))*drdT*w(i)))

 30        CONTINUE

      dPPTT = dP0dTT + dPdTT + P11*ddelMT - D/Dc*dPdMTT/d2PdM2 +
     1        (P11 + d2PdMT)*(D/Dc*alpha/d2PdM2 + 
     2        D/Dc*dPdMMT/d2PdM2**2)

      pterm = P/Pc + dPdTcd

*** compute (d(alpha)/dT)P

      dalLVS  = Tc*Dc**2/D**2/T**2 * (-2.0d0/T*d2PdM2*pterm + 
     1          2.0d0*alpha*d2PdM2*pterm + pterm*dPdMMT +
     2          d2PdM2*dPPTT)

      RETURN

      END

*********************************************************************

*** backup - Save Pfind COMMON values during saturation check.

      SUBROUTINE backup

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION s(2), xk(2)

      COMMON /satur/ Dliq, Dvap, DH2O, iphase
      COMMON /param/ r1, th1
      COMMON /deriv/ amu, s, Pw, Tw, dTw, dM0dT, dP0dT,
     1               d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
      COMMON /deri2/ dPdD, dPdT
      COMMON /store/ sav2, sav3, sav4, sav5, sav6, sav7, sav8, 
     1               sav9, sav10, sav11, sav12, sav13, sav14, sav15, 
     2               sav16, sav17, sav18, sav19, isav1

      SAVE


      isav1 = iphase
      
      sav2  = r1
      sav3  = th1

      sav4  = amu
      sav5  = s(1)
      sav6  = s(2)
      sav7  = Pw
      sav8  = Tw
      sav9  = dTw
      sav10 = dM0dT
      sav11 = dP0dT
      sav12 = d2PdM2
      sav13 = d2PdMT
      sav14 = d2PdT2
      sav15 = p0th
      sav16 = p1th
      sav17 = xk(1)
      sav18 = xk(2)

      sav19 = dPdD

      RETURN
      END

*********************************************************************

*** restor - Restore Pfind COMMON values after saturation check.


      SUBROUTINE restor

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION s(2), xk(2)

      COMMON /satur/ Dliq, Dvap, DH2O, iphase
      COMMON /param/ r1, th1
      COMMON /deriv/ amu, s, Pw, Tw, dTw, dM0dT, dP0dT,
     1               d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
      COMMON /deri2/ dPdD, dPdT
      COMMON /store/ sav2, sav3, sav4, sav5, sav6, sav7, sav8, 
     1               sav9, sav10, sav11, sav12, sav13, sav14, sav15, 
     2               sav16, sav17, sav18, sav19, isav1

      SAVE

       
      iphase = isav1

      r1     = sav2
      th1    = sav3
    
      amu    = sav4
      s(1)   = sav5
      s(2)   = sav6
      Pw     = sav7
      Tw     = sav8
      dTw    = sav9
      dM0dT  = sav10
      dP0dT  = sav11
      d2PdM2 = sav12
      d2PdMT = sav13
      d2PdT2 = sav14
      p0th   = sav15
      p1th   = sav16
      xk(1)  = sav17
      xk(2)  = sav18

      dPdD   = sav19

      RETURN
      END

**********************************************************************

*** load - Load thermodynamic and transport property values from
*          ptemp into props.

      SUBROUTINE load(phase,ptemp,props)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP  = 23, NPROP2 = 46)

      DOUBLE PRECISION  ptemp(NPROP), props(NPROP2)
      INTEGER           phase, key(NPROP,2)

      SAVE

      DATA  key
     1   /  1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23,
     2     25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 
     3      2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 
     4     26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46  /


      DO 10  i = 1,NPROP 
 10        props(key(i,phase)) = ptemp(i)

      RETURN
      END

******************************************************************

*** tpset - Dimension triple point  U, S, H, A, G  values (in J/g from
*           Table 2, Helgeson & Kirkham, 1974a) into user-specified units.

      SUBROUTINE tpset

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
      COMMON /tpoint/ Utripl, Stripl, Htripl, Atripl, Gtripl, 
     1                Ttripl, Ptripl, Dltrip, Dvtrip

      SAVE

      DATA       Utr,        Str,       Htr,        Atr,        Gtr
     1     / -15766.0d0,  3.5144d0, -15971.0d0, -12870.0d0, -13073.0d0 /

      
      Utripl = Utr * fh
      Stripl = Str * fh
      Htripl = Htr * fh
      Atripl = Atr * fh
      Gtripl = Gtr * fh

      END

****************************************************************************

*** triple - Convert  U, S, H, A, G  values computed with reference to  
*            zero triple point properties (Haar et al., 1984; 
*            Levelt Sengers et al., 1983) into values referenced to 
*            triple point properties given by  Helgeson and Kirkham, 1974a.

      SUBROUTINE triple(T,wpzero)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP = 23)

      DOUBLE PRECISION  wpzero(NPROP)
      INTEGER  A, G, S, U, H

      COMMON /tpoint/ Utr, Str, Htr, Atr, Gtr, 
     1                Ttripl, Ptripl, Dltrip, Dvtrip

      SAVE

      DATA    A, G, S, U, H
     1      / 1, 2, 3, 4, 5 /


      wpzero(S) = wpzero(S) + Str

      TS = T*wpzero(S) - Ttripl*Str

      wpzero(G) = wpzero(H) - TS + Gtr            
      wpzero(A) = wpzero(U) - TS + Atr 

      wpzero(H) = wpzero(H) + Htr
      wpzero(U) = wpzero(U) + Utr

      END

*********************************************************************

*** power - Returns  base**exp  utilizing the intrinsic FORTRAN 
*           exponentiation function in such a manner so as to 
*           insure computation of machine-independent values 
*           for all defined exponentiations.  Attempted undefined
*           exponentiations produce an error message and cause
*           program termination.

      DOUBLE PRECISION FUNCTION power(base,exp)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER rterm, wterm, reacf, pronf, tabf, plotf(6)

      COMMON /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf

      SAVE

      DATA TOL / 1.0d-7 /

      
      IF (base .GT. 0.0d0) THEN
           power = base**exp
      ELSE
           IF (DABS(base) .GT. TOL) THEN
                IF (DBLE(INT(exp)) .NE. exp) THEN
                     WRITE(wterm,10) base, exp
 10                  FORMAT(/,' neg base ** real exp is complex',
     1                      /,' base,exp: ',2e20.13,/)
                     STOP
                ELSE
                     IF (MOD(exp,2.0d0) .EQ. 0.0d0) THEN
                          power =  (-base)**exp
                     ELSE
                          power = -((-base)**exp)
                     END IF
                END IF
           ELSE
                IF (exp .GT. 0.0d0) THEN
                     power = 0.0d0
                ELSE
                     WRITE(wterm,20) base, exp
 20                  FORMAT(/,' zero base ** (exp <= 0) is undefined',
     1                      /,' base,exp: ',2e20.13)
                     STOP
                END IF
           END IF
      END IF

      RETURN
      END   

***********************************************************************

*** TdegK - Returns input temperature  t  converted from 
*           user-specified units to degK.

      DOUBLE PRECISION FUNCTION TdegK(it,t)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      SAVE


      GO TO (1,2,3,4), it

 1    TdegK = t
      RETURN

 2    TdegK = t + 273.15d0
      RETURN

 3    TdegK = t / 1.8d0
      RETURN

 4    TdegK = (t + 459.67d0) / 1.8d0
      RETURN

      END

***********************************************************************

*** TdegUS - Returns input temperature  t  converted 
*            from degK to user-specified units.

      DOUBLE PRECISION FUNCTION TdegUS(it,t)
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      SAVE


      GO TO (1,2,3,4), it

 1    TdegUS = t
      RETURN

 2    TdegUS = t - 273.15d0
      RETURN

 3    TdegUS = t * 1.8d0
      RETURN

 4    TdegUS = t * 1.8d0 - 459.67d0
      RETURN

      END

*********************************************************************

*** dim[HGK,LVS] - Dimensioning routines for H2O88.

*********************************************************************

*** dimHGK - Dimensions thermodynamic and transport property values 
*            computed from the HGK equation of state per user-specified 
*            choice of units.     
            
      SUBROUTINE dimHGK(isat,itripl,t,p,d,epseqn)
                                  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP = 23)

      DOUBLE PRECISION  wprops(NPROP), wpliq(NPROP)
      INTEGER           aw, gw, sw, uw, hw, cvw, cpw, vsw, alw, bew,
     1                  diw, viw, tcw, stw, tdw, Prw, vikw, albew,
     2                  ZBw, YBw, QBw, dalwdT, XBw
      INTEGER  epseqn    

      COMMON /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
      COMMON /fcts/   ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, dpdd,
     1                cjtt, cjth
      COMMON /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
      COMMON /RTcurr/ rt    
      COMMON /wpvals/ wprops, wpliq

      SAVE

      DATA aw, gw, sw, uw, hw, cvw, cpw, vsw, alw, bew, diw, viw,
     1     tcw, stw, tdw, Prw, vikw, albew, ZBw, YBw, QBw, 
     2     dalwdT, XBw
     3   /  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12,
     4     13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 /


      wprops(aw)   = ad * rt * fh
      wprops(gw)   = gd * rt * fh
      wprops(sw)   = sd * gascon * fh * ft
      wprops(uw)   = ud * rt * fh
      wprops(hw)   = hd * rt * fh
      wprops(cvw)  = cvd * gascon * fh * ft
      wprops(cpw)  = cpd * gascon * fh * ft
      wprops(vsw)  = DSQRT(DABS(cpd*dpdd*1.0d3/cvd)) * fs
      wprops(bew)  = 1.0d0 / (d * dpdd * fp)
      wprops(alw)  = d * dvdt
      wprops(dalwdT) = dalHGK(d,t,wprops(alw))

 
      pbars = p*1.0d1
      dkgm3 = d * 1.0d3   
      betaPa = wprops(bew)*fp / 1.0d6
      betab  = wprops(bew)*fp / 1.0d1
      CpJKkg = wprops(cpw)/fh/ft * 1.0d3

      wprops(viw)  = viscos(t,pbars,dkgm3,betaPa) * fvd
      wprops(tcw)  = thcond(t,pbars,dkgm3,wprops(alw),betaPa) * fc * ft
      IF ((isat .EQ. 0) .OR. (isat .EQ. 2)) THEN
           wprops(stw) = 0.0d0
      ELSE
           wprops(stw) = surten(t) * fst
      END IF

      CALL Born92(t,pbars,dkgm3/1.0d3,betab,wprops(alw),wprops(dalwdT),
     1            wprops(diw),wprops(ZBw),wprops(QBw),wprops(YBw),
     2            wprops(XBw),epseqn)

      wprops(tdw)   = wprops(tcw)/fc/ft  / (dkgm3 * CpJKkg) * fvk
      IF (wprops(tcw) .NE. 0.0d0) THEN
           wprops(Prw) = wprops(viw)/fvd * CpJKkg / (wprops(tcw)/fc/ft)
      ELSE
           wprops(Prw) = 0.0d0
      END IF
      wprops(vikw)  = wprops(viw)/fvd / dkgm3 * fvk
      wprops(albew) = wprops(alw) / wprops(bew)

      IF (itripl .EQ. 1) CALL triple(t,wprops)

      END 

*****************************************************************************

*** dimLVS - Dimension critical region properties per user-specs
*            and load into tprops.

      SUBROUTINE dimLVS(isat,itripl,theta,T,Pbars,dl,dv,tprops,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NPROP = 23)

      DOUBLE PRECISION  tprops(NPROP)
      INTEGER   aw, gw, sw, uw, hw, cvw, cpw, vsw, alw, bew,
     1          diw, viw, tcw, stw, tdw, Prw, vikw, albew,
     2          ZBw, YBw, QBw, dalwdT, XBw
      INTEGER  epseqn


      COMMON /therm/   AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw,
     1                 heat, Speed
      COMMON /satur/   Dliq, Dvap, DH2O, iphase
      COMMON /units/   ft, fd, fvd, fvk, fs, fp, fh, fst, fc
*****************************************************************
      COMMON /abc2/    r, th
*****************************************************************

      SAVE

      DATA aw, gw, sw, uw, hw, cvw, cpw, vsw, alw, bew, diw, viw,
     1     tcw, stw, tdw, Prw, vikw, albew, ZBw, YBw, QBw,
     2     dalwdT, XBw
     3   /  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  
     4     13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 /

  
      IF (isat .EQ. 1) THEN
          dv   = Dvap
          dl   = Dliq
      END IF

      tprops(aw)  = AE * fh
      tprops(gw)  = GE * fh
      tprops(sw)  = Entrop * fh * ft
      tprops(uw)  = U * fh
      tprops(hw)  = H * fh
      tprops(cvw) = Cv * fh * ft
      tprops(cpw) = Cp * fh * ft
      tprops(vsw) = Speed * fs
      tprops(bew) = betaw / fp
      tprops(alw) = alphw
*****************************************************************
      th = theta
      tprops(dalwdT) = dalLVS(DH2O,T,Pbars/1.0d1,tprops(alw))
*****************************************************************

      CpJKkg  = Cp * 1.0d3
      betaPa  = betaw / 1.0d6
      betab   = betaw / 1.0d1

      IF (DABS(theta) .NE. 1.0d0) THEN
           dkgm3 = DH2O
           tprops(stw) = 0.0d0
      ELSE
           IF (theta .LT. 0.0d0) THEN
                dkgm3 = Dvap
                tprops(stw) = 0.0d0
           ELSE
                dkgm3 = Dliq
                dkgm3 = Dliq
                tprops(stw) = surten(T) * fst
           END IF
      END IF

      CALL Born92(T,Pbars,dkgm3/1.0d3,betab,tprops(alw),tprops(dalwdT),
     1            tprops(diw),tprops(ZBw),tprops(QBw),tprops(YBw),
     2            tprops(XBw),epseqn)

      tprops(viw)  = viscos(T,Pbars,dkgm3,betaPa) * fvd
      tprops(tcw)  = thcond(T,Pbars,dkgm3,tprops(alw),betaPa) * fc * ft

      tprops(tdw)  = tprops(tcw)/fc/ft  / (dkgm3 * CpJKkg) * fvk
      tprops(Prw)  = tprops(viw)/fvd * CpJKkg / (tprops(tcw)/fc/ft)
      tprops(vikw) = tprops(viw)/fvd / dkgm3 * fvk
      tprops(albew) = tprops(alw) / tprops(bew)

      IF (itripl .EQ. 1)  CALL triple(T,tprops)

      END      

**********************************************************************

*** tran88 - Set of FORTRAN77 functions that compute transport
*            properties of fluid H2O.  Input state parameters 
*            should be computed from the Haar et al. (1984)
*            and Levelt Sengers et al. (1983) equations of state in 
*            order to facilitate comparision with published tabular
*            values referenced below for each function.
*
**********************************************************************

***   programmer:  James W. Johnson
***   abandoned:   20 January 1988

**********************************************************************

*** viscos - Returns dynamic viscosity of H2O in kg/m*s (= Pa*s)
*            if  Tk, Pbars  falls within the validity region (specified 
*            by the initial IF statement) of the Watson et al. (1980)
*            equation; otherwise returns zero.  See equations 3.1-2 and
*            4.1-5 and Tables 1, 6, and 8 from Sengers and 
*            Kamgar-Parsi (1984).

      DOUBLE PRECISION FUNCTION viscos(Tk,Pbars,Dkgm3,betaPa)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (Tstar = 647.270d0)
      PARAMETER (Dstar = 317.763d0)
      PARAMETER (Pstar = 22.1150d6)
      PARAMETER (ustar = 1.0d-6)

      DOUBLE PRECISION  a(4), b(6,7)

      SAVE

      DATA a / 0.0181583d0,  0.0177624d0,  0.0105287d0, -0.0036744d0 /

      DATA b / 0.5132047d0,  0.3205656d0,  0.0d0,        0.0d0,
     1        -0.7782567d0,  0.1885447d0,  0.2151778d0,  0.7317883d0,
     2         1.2410440d0,  1.4767830d0,  0.0d0,        0.0d0,
     3        -0.2818107d0, -1.0707860d0, -1.2631840d0,  0.0d0,
     4         0.0d0,        0.0d0,        0.1778064d0,  0.4605040d0,
     5         0.2340379d0, -0.4924179d0,  0.0d0,        0.0d0,
     6        -0.0417661d0,  0.0d0,        0.0d0,        0.1600435d0,
     7         0.0d0,        0.0d0,        0.0d0,       -0.01578386d0,
     8         0.0d0,        0.0d0,        0.0d0,        0.0d0,
     9         0.0d0,        0.0d0,        0.0d0,      -0.003629481d0,
     1         0.0d0,        0.0d0 /
 
      DATA TOL /1.0d-2/


      viscos = 0.0d0
      TdegC  = Tk - 273.15d0

      IF ((Pbars .GT. 5000.0d0+TOL) .OR.
     1   ((Pbars .GT. 3500.0d0+TOL).AND.(TdegC .GT. 150.0d0+TOL)).OR. 
     2   ((Pbars .GT. 3000.0d0+TOL).AND.(TdegC .GT. 600.0d0+TOL)) .OR.
     3   (TdegC  .GT. 900.0d0+TOL))  RETURN

      T = Tk / Tstar
      D = Dkgm3 / Dstar
      
      sum = 0.0d0
      DO 10  i=0,3
 10        sum = sum + a(i+1)/T**i
      u0 = ustar * DSQRT(T) / sum

      sum = 0.0d0
      DO 20  i=0,5
           DO 20  j=0,6
 20             sum = sum + b(i+1,j+1) * (1.0d0/T-1)**i * (D-1)**j
      u1 = DEXP(D*sum)

      IF ((0.997d0 .LE. T) .AND. (T .LE. 1.0082d0) .AND.
     1    (0.755d0 .LE. D) .AND. (D .LE. 1.2900d0)) THEN
           xt = Pstar/Dstar**2 * betaPa * Dkgm3**2
           IF (xt .LT. 22.0d0) THEN
                u2 = 1.0d0
           ELSE
                u2 = 0.922 * power(xt,0.0263d0)
           END IF
      ELSE
           u2 = 1.0d0
      END IF

      viscos = u0 * u1 * u2

      RETURN
      END

*****************************************************************

*** thcond - Returns thermal conductivity of H2O in J/m*deg*s (=W/m*deg)
*            if  Tk, Pbars  falls within the validity region (specified 
*            by the initial IF statement) of the Sengers et al. (1984)
*            equation; returns zero otherwise.  See equations 3.2-14
*            and tables 2-5 and I.5-6 from the above reference.

      DOUBLE PRECISION FUNCTION thcond(Tk,Pbars,Dkgm3,alph,betaPa)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (Tstar = 647.270d0)
      PARAMETER (Dstar = 317.763d0)
      PARAMETER (Pstar = 22.1150d6)
      PARAMETER (ustar = 1.0d-6)
      PARAMETER (C     = 3.7711d-8)

      DOUBLE PRECISION  aL(4), au(4), bL(6,5), bu(5,6), L0, L1, L2

      SAVE

      DATA aL / 0.2022230d1,  0.1411166d2,  0.5255970d1, -0.2018700d1 /

      DATA au / 0.0181583d0,  0.0177624d0,  0.0105287d0, -0.0036744d0 /

      DATA bL / 1.329304600d0, -0.404524370d0,  0.244094900d0,
     1          0.018660751d0, -0.129610680d0,  0.044809953d0,
     2          1.701836300d0, -2.215684500d0,  1.651105700d0,
     3         -0.767360020d0,  0.372833440d0, -0.112031600d0,
     4          5.224615800d0, -1.012411100d1,  4.987468700d0,
     5         -0.272976940d0, -0.430833930d0,  0.133338490d0,
     6          8.712767500d0, -9.500061100d0,  4.378660600d0,
     7         -0.917837820d0,  0.0d0,          0.0d0,
     8         -1.852599900d0,  0.934046900d0,  0.0d0,
     9          0.0d0,          0.0d0,          0.0d0  /

      DATA bu / 0.5019380d0,  0.2356220d0, -0.2746370d0,  0.1458310d0,
     1         -0.0270448d0,  0.1628880d0,  0.7893930d0, -0.7435390d0,
     2          0.2631290d0, -0.0253093d0, -0.1303560d0,  0.6736650d0,
     3         -0.9594560d0,  0.3472470d0, -0.0267758d0,  0.9079190d0,
     4          1.2075520d0, -0.6873430d0,  0.2134860d0, -0.0822904d0,
     5         -0.5511190d0,  0.0670665d0, -0.4970890d0,  0.1007540d0,
     6          0.0602253d0,  0.1465430d0, -0.0843370d0,  0.1952860d0,
     7         -0.0329320d0, -0.0202595d0  /

      DATA TOL /1.0d-2/


      thcond = 0.0d0
      TdegC  = Tk - 273.15d0

      IF ((Pbars .GT. 4000.0d0+TOL) .OR.
     1   ((Pbars .GT. 2000.0d0+TOL).AND.(TdegC .GT. 125.0d0+TOL)).OR.
     2   ((Pbars .GT. 1500.0d0+TOL).AND.(TdegC .GT. 400.0d0+TOL)).OR.
     3   (TdegC  .GT. 800.0d0+TOL))  RETURN

      T = Tk / Tstar
      D = Dkgm3 / Dstar

      sum = 0.0d0
      DO 10  i=0,3
 10        sum = sum + aL(i+1)/T**i
      L0 = DSQRT(T) / sum

      sum = 0.0d0
      DO 20  i=0,4
           DO 20  j=0,5
 20             sum = sum + bL(j+1,i+1) * (1.0d0/T-1)**i * (D-1)**j
      L1 = DEXP(D*sum)

      sum = 0.0d0
      DO 40  i=0,3
 40        sum = sum + au(i+1)/T**i
      u0 = ustar * DSQRT(T) / sum

      sum = 0.0d0
      DO 50  i=0,5
           DO 50  j=0,4
 50             sum = sum + bu(j+1,i+1) * (1.0d0/T-1)**i * (D-1)**j
      u1 = DEXP(D*sum)

      xt   = Pstar/Dstar**2 * betaPa * Dkgm3**2
      dPdT = Tstar/Pstar * alph/betaPa

      L2 = C / (u0*u1) * (T/D)**2 * dPdT**2 * power(xt,0.4678d0) * 
     1     DSQRT(D) * DEXP(-18.66d0*(T-1)**2 - (D-1)**4)

      thcond = L0 * L1 + L2

      RETURN
      END

******************************************************************

*** surten - Returns the surface tension of vapor-saturated liquid
*            H2O in MPa*cm (converted from N/m) as computed from 
*            the Vargaftik et al. (1983) equation.  See equations
*            10.1-2, Kestin et al. (1984); compare also equation
*            C.5 and table 11, Haar et al. (1984).

      DOUBLE PRECISION FUNCTION surten(Tsatur)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (Ttripl = 273.16d0)
      PARAMETER (Tcrit  = 647.067d0)
      PARAMETER (Tstar  = 647.27d0)
      PARAMETER (Tcstar = 0.999686d0)
      PARAMETER (v =  1.256d0)
      PARAMETER (B = -0.625d0)
      PARAMETER (stref = 0.2358d0)
      PARAMETER (FPTOL = 1.0d-10)

      SAVE


      IF ((Tsatur .LT. Ttripl) .OR. (Tsatur .GT. Tcrit)) THEN
           surten = 0.0d0
           RETURN
      END IF

      IF (Tsatur .GE. Tcrit-FPTOL) THEN
           Tnorm = 0.0d0
      ELSE
           Tnorm = (Tcstar - Tsatur/Tstar) / Tcstar
      END IF

      surten = stref * power(Tnorm,v) * (1.0d0 + B*Tnorm)

      RETURN

      END

******************************************************************

*** Born92 - Computes the Z, Q, Y, and X Born functions at TK, Pbars.
***        
***        epseqn = 1 ...... use Helgeson-Kirkham (1974) equation 
***        epseqn = 2 ...... use Pitzer (1983) equation
***        epseqn = 3 ...... use Uematsu-Franck (1980) equation
***        epseqn = 4 ...... use Johnson-Norton (1991) equation
***        epseqn = 5 ...... use Archer-Wang (1990) equation
***                         
      SUBROUTINE Born92(TK,Pbars,Dgcm3,betab,alphaK,daldT,
     1                  eps,Z,Q,Y,X,epseqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (TMAX = 1000.0d0, PMAX = 5000.0d0, TOL = 1.0d-3)

      INTEGER  epseqn

      SAVE


      eps = 0.0d0
      Z   = 0.0d0
      Y   = 0.0d0
      Q   = 0.0d0
      X   = 0.0d0

      TdegC = TK - 273.15d0

***   The following line can be commented out to facilitate probably
***   unreliable, yet potentially useful, predictive calculations 
***   at state conditions beyond the validity limits of the aqueous 
***   species equation of state.

      IF ((TdegC .GT. TMAX+TOL) .OR. (Pbars .GT. PMAX+TOL)) RETURN

*      IF (epseqn .EQ. 1) THEN
*           CALL HK74(TK,Dgcm3,betab,alphaK,daldT,
*     1               eps,dedP,dedT,d2edT2)
*           CALL epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
*           RETURN
*      END IF

*      IF (epseqn .EQ. 2) THEN
*           CALL Pitz83(TK,Dgcm3,betab,alphaK,daldT,
*     1                 eps,dedP,dedT,d2edT2)
*           CALL epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
*           RETURN
*      END IF

*      IF (epseqn .EQ. 3) THEN
*           CALL UF80(TK,Dgcm3,betab,alphaK,daldT,
*     1               eps,dedP,dedT,d2edT2)
*           CALL epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
*           RETURN
*      END IF
     
      IF (epseqn .EQ. 4) THEN
           CALL JN91(TK,Dgcm3,betab,alphaK,daldT,
     1                eps,dedP,dedT,d2edT2)
           CALL epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
           RETURN
      END IF

*      IF (epseqn .EQ. 5) THEN
*           Dkgm3 = Dgcm3 * 1.0d3
*           PMPa  = Pbars / 1.0d1
*           betam = betab * 1.0d1
*           CALL AW90(TK,PMPa,Dkgm3,betam,alphaK,daldT,
*     1                 eps,dedP,dedT,d2edT2)
****        convert  dedP  FROM  MPa**-1  TO  bars**-1
*           dedP = dedP / 1.0d1
*           CALL epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
*           RETURN
*      END IF
     
      END

*********************************************************************

*** JN91 - Compute (eps, dedP, dedT, d2edT2)(T,D) using equations 
***        given by Johnson and Norton (1991); fit parameters 
***        regressed from least squares fit to dielectric data
***        consistent with the HK74 equation and low temperatures,
***        and with the Pitz83 equation at high temperatures.
***
***          Units: T ............... K
***                 D ............... g/cm**3
***                 beta, dedP ...... bar**(-1)
***                 alpha, dedT ..... K**(-1)
***                 daldT, d2edT2 ... K**(-2)


      SUBROUTINE JN91(T,D,beta,alpha,daldT,eps,dedP,dedT,d2edT2)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION  a(10), c(5), dcdT(5), dc2dTT(5)

      SAVE

      DATA Tref / 298.15d0 /

      DATA a / 
     1          0.1470333593E+02, 
     2          0.2128462733E+03, 
     3         -0.1154445173E+03, 
     4          0.1955210915E+02, 
     5         -0.8330347980E+02, 
     6          0.3213240048E+02, 
     7         -0.6694098645E+01, 
     8         -0.3786202045E+02, 
     9          0.6887359646E+02, 
     1         -0.2729401652E+02 /

      Tn = T / Tref

      c(1)      = 1.0d0
      dcdT(1)   = 0.0d0
      dc2dTT(1) = 0.0d0

      c(2)      = a(1)/Tn
      dcdT(2)   = -a(1)*Tref/T**2
      dc2dTT(2) = 2.0d0*a(1)*Tref/T**3

      c(3)      = a(2)/Tn + a(3) + a(4)*Tn           
      dcdT(3)   = -a(2)*Tref/T**2 + a(4)/Tref           
      dc2dTT(3) = 2.0d0*a(2)*Tref/T**3

      c(4)      = a(5)/Tn + a(6)*Tn + a(7)*Tn**2
      dcdT(4)   = -a(5)*Tref/T**2 + a(6)/Tref 
     1            + 2.0d0*a(7)*T/Tref**2
      dc2dTT(4) = 2.0d0*a(5)*Tref/T**3 + 2.0d0*a(7)/Tref**2

      c(5)      = a(8)/Tn**2 + a(9)/Tn + a(10)
      dcdT(5)   = -2.0d0*a(8)*Tref**2/T**3 - a(9)*Tref/T**2
      dc2dTT(5) = 6.0d0*a(8)*Tref**2/T**4 + 2.0d0*a(9)*Tref/T**3
    
      eps = 0.0d0
      DO 50 k=1,5
   50      eps = eps + c(k)*D**(k-1)

      dedP = 0.0d0
      DO 100  j = 0,4
  100      dedP = dedP + j*c(j+1)*D**j
      dedP = beta * dedP

      dedT = 0.0d0
      DO 200  j = 0,4
  200      dedT = dedT + D**j*(dcdT(j+1) - j*alpha*c(j+1))

      d2edT2 = 0.0d0
      DO 300  j = 0,4
  300      d2edT2 = d2edT2 + D**j*(dc2dTT(j+1) - j*(alpha*dcdT(j+1) +
     1         c(j+1)*daldT) - j*alpha*(dcdT(j+1) - j*alpha*c(j+1)))

      END

***************************************************************

*** epsBrn - Compute the Z, Q, Y, and X Born functions from their
***          eps, dedP, dedT, and d2edT2 counterparts.

      SUBROUTINE epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      SAVE

      Z = -1.0d0/eps
      Q =  1.0d0/eps**2 * dedP
      Y =  1.0d0/eps**2 * dedT
      X =  1.0d0/eps**2 * d2edT2 - 2.0d0*eps*Y**2
       
      END
