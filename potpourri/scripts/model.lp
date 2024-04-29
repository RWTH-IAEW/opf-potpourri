\* Source Pyomo model name=unknown *\

min 
OBJ:
-9000000 alpha(D10)
-3500000.0000000005 alpha(D11)
-6100000 alpha(D12)
-13500000 alpha(D13)
-14899999.999999998 alpha(D14)
-16940000 alpha(D2)
-94200000 alpha(D3)
-47800000 alpha(D4)
-7600000 alpha(D5)
-11200000 alpha(D6)
-29500000 alpha(D9)
+100 pDis(B1_0)
+100 pDis(B1_1)
+100 pDis(B1_2)
+100 pDis(B1_3)
+100 pDis(B1_4)
+100 pDis(B1_5)
+100 pDis(B1_6)
+100 pDis(B1_7)
+100 pDis(B1_8)
+100 pDis(B1_9)
+100 pDis(B2_0)
+100 pDis(B2_1)
+100 pDis(B2_2)
+100 pDis(B2_3)
+100 pDis(B2_4)
+100 pDis(B2_5)
+100 pDis(B2_6)
+100 pDis(B2_7)
+100 pDis(B2_8)
+100 pDis(B2_9)
+2000 pG(G1_0)
+2000 pG(G1_1)
+2000 pG(G1_2)
+2000 pG(G1_3)
+2000 pG(G1_4)
+2000 pG(G1_5)
+2000 pG(G1_6)
+2000 pG(G1_7)
+2000 pG(G1_8)
+2000 pG(G1_9)
+2000 pG(G2_0)
+2000 pG(G2_1)
+2000 pG(G2_2)
+2000 pG(G2_3)
+2000 pG(G2_4)
+2000 pG(G2_5)
+2000 pG(G2_6)
+2000 pG(G2_7)
+2000 pG(G2_8)
+2000 pG(G2_9)
+4000 pG(G3_0)
+4000 pG(G3_1)
+4000 pG(G3_2)
+4000 pG(G3_3)
+4000 pG(G3_4)
+4000 pG(G3_5)
+4000 pG(G3_6)
+4000 pG(G3_7)
+4000 pG(G3_8)
+4000 pG(G3_9)
+4000 pG(G4_0)
+4000 pG(G4_1)
+4000 pG(G4_2)
+4000 pG(G4_3)
+4000 pG(G4_4)
+4000 pG(G4_5)
+4000 pG(G4_6)
+4000 pG(G4_7)
+4000 pG(G4_8)
+4000 pG(G4_9)
+4000 pG(G5_0)
+4000 pG(G5_1)
+4000 pG(G5_2)
+4000 pG(G5_3)
+4000 pG(G5_4)
+4000 pG(G5_5)
+4000 pG(G5_6)
+4000 pG(G5_7)
+4000 pG(G5_8)
+4000 pG(G5_9)
+ [
+860.58519799999988 pG(G1_0) ^ 2
+860.58519799999988 pG(G1_1) ^ 2
+860.58519799999988 pG(G1_2) ^ 2
+860.58519799999988 pG(G1_3) ^ 2
+860.58519799999988 pG(G1_4) ^ 2
+860.58519799999988 pG(G1_5) ^ 2
+860.58519799999988 pG(G1_6) ^ 2
+860.58519799999988 pG(G1_7) ^ 2
+860.58519799999988 pG(G1_8) ^ 2
+860.58519799999988 pG(G1_9) ^ 2
+5000 pG(G2_0) ^ 2
+5000 pG(G2_1) ^ 2
+5000 pG(G2_2) ^ 2
+5000 pG(G2_3) ^ 2
+5000 pG(G2_4) ^ 2
+5000 pG(G2_5) ^ 2
+5000 pG(G2_6) ^ 2
+5000 pG(G2_7) ^ 2
+5000 pG(G2_8) ^ 2
+5000 pG(G2_9) ^ 2
+200 pG(G3_0) ^ 2
+200 pG(G3_1) ^ 2
+200 pG(G3_2) ^ 2
+200 pG(G3_3) ^ 2
+200 pG(G3_4) ^ 2
+200 pG(G3_5) ^ 2
+200 pG(G3_6) ^ 2
+200 pG(G3_7) ^ 2
+200 pG(G3_8) ^ 2
+200 pG(G3_9) ^ 2
+200 pG(G4_0) ^ 2
+200 pG(G4_1) ^ 2
+200 pG(G4_2) ^ 2
+200 pG(G4_3) ^ 2
+200 pG(G4_4) ^ 2
+200 pG(G4_5) ^ 2
+200 pG(G4_6) ^ 2
+200 pG(G4_7) ^ 2
+200 pG(G4_8) ^ 2
+200 pG(G4_9) ^ 2
+200 pG(G5_0) ^ 2
+200 pG(G5_1) ^ 2
+200 pG(G5_2) ^ 2
+200 pG(G5_3) ^ 2
+200 pG(G5_4) ^ 2
+200 pG(G5_5) ^ 2
+200 pG(G5_6) ^ 2
+200 pG(G5_7) ^ 2
+200 pG(G5_8) ^ 2
+200 pG(G5_9) ^ 2
] / 2
+254240000 ONE_VAR_CONSTANT

s.t.

c_e_KCL_real(1_0)_:
-1 pChar(B1_0)
-1 pD(D1_0)
+1 pDis(B1_0)
+1 pG(G1_0)
-1 pLfrom(L1_12_0)
-1 pLfrom(L2_15_0)
= 0

c_e_KCL_real(1_1)_:
-1 pChar(B1_1)
-1 pD(D1_1)
+1 pDis(B1_1)
+1 pG(G1_1)
-1 pLfrom(L1_12_1)
-1 pLfrom(L2_15_1)
= 0

c_e_KCL_real(1_2)_:
-1 pChar(B1_2)
-1 pD(D1_2)
+1 pDis(B1_2)
+1 pG(G1_2)
-1 pLfrom(L1_12_2)
-1 pLfrom(L2_15_2)
= 0

c_e_KCL_real(1_3)_:
-1 pChar(B1_3)
-1 pD(D1_3)
+1 pDis(B1_3)
+1 pG(G1_3)
-1 pLfrom(L1_12_3)
-1 pLfrom(L2_15_3)
= 0

c_e_KCL_real(1_4)_:
-1 pChar(B1_4)
-1 pD(D1_4)
+1 pDis(B1_4)
+1 pG(G1_4)
-1 pLfrom(L1_12_4)
-1 pLfrom(L2_15_4)
= 0

c_e_KCL_real(1_5)_:
-1 pChar(B1_5)
-1 pD(D1_5)
+1 pDis(B1_5)
+1 pG(G1_5)
-1 pLfrom(L1_12_5)
-1 pLfrom(L2_15_5)
= 0

c_e_KCL_real(1_6)_:
-1 pChar(B1_6)
-1 pD(D1_6)
+1 pDis(B1_6)
+1 pG(G1_6)
-1 pLfrom(L1_12_6)
-1 pLfrom(L2_15_6)
= 0

c_e_KCL_real(1_7)_:
-1 pChar(B1_7)
-1 pD(D1_7)
+1 pDis(B1_7)
+1 pG(G1_7)
-1 pLfrom(L1_12_7)
-1 pLfrom(L2_15_7)
= 0

c_e_KCL_real(1_8)_:
-1 pChar(B1_8)
-1 pD(D1_8)
+1 pDis(B1_8)
+1 pG(G1_8)
-1 pLfrom(L1_12_8)
-1 pLfrom(L2_15_8)
= 0

c_e_KCL_real(1_9)_:
-1 pChar(B1_9)
-1 pD(D1_9)
+1 pDis(B1_9)
+1 pG(G1_9)
-1 pLfrom(L1_12_9)
-1 pLfrom(L2_15_9)
= 0

c_e_KCL_real(2_0)_:
-1 pChar(B2_0)
-1 pD(D2_0)
+1 pDis(B2_0)
+1 pG(G2_0)
-1 pLfrom(L3_23_0)
-1 pLfrom(L4_24_0)
-1 pLfrom(L5_25_0)
-1 pLto(L1_12_0)
= 0

c_e_KCL_real(2_1)_:
-1 pChar(B2_1)
-1 pD(D2_1)
+1 pDis(B2_1)
+1 pG(G2_1)
-1 pLfrom(L3_23_1)
-1 pLfrom(L4_24_1)
-1 pLfrom(L5_25_1)
-1 pLto(L1_12_1)
= 0

c_e_KCL_real(2_2)_:
-1 pChar(B2_2)
-1 pD(D2_2)
+1 pDis(B2_2)
+1 pG(G2_2)
-1 pLfrom(L3_23_2)
-1 pLfrom(L4_24_2)
-1 pLfrom(L5_25_2)
-1 pLto(L1_12_2)
= 0

c_e_KCL_real(2_3)_:
-1 pChar(B2_3)
-1 pD(D2_3)
+1 pDis(B2_3)
+1 pG(G2_3)
-1 pLfrom(L3_23_3)
-1 pLfrom(L4_24_3)
-1 pLfrom(L5_25_3)
-1 pLto(L1_12_3)
= 0

c_e_KCL_real(2_4)_:
-1 pChar(B2_4)
-1 pD(D2_4)
+1 pDis(B2_4)
+1 pG(G2_4)
-1 pLfrom(L3_23_4)
-1 pLfrom(L4_24_4)
-1 pLfrom(L5_25_4)
-1 pLto(L1_12_4)
= 0

c_e_KCL_real(2_5)_:
-1 pChar(B2_5)
-1 pD(D2_5)
+1 pDis(B2_5)
+1 pG(G2_5)
-1 pLfrom(L3_23_5)
-1 pLfrom(L4_24_5)
-1 pLfrom(L5_25_5)
-1 pLto(L1_12_5)
= 0

c_e_KCL_real(2_6)_:
-1 pChar(B2_6)
-1 pD(D2_6)
+1 pDis(B2_6)
+1 pG(G2_6)
-1 pLfrom(L3_23_6)
-1 pLfrom(L4_24_6)
-1 pLfrom(L5_25_6)
-1 pLto(L1_12_6)
= 0

c_e_KCL_real(2_7)_:
-1 pChar(B2_7)
-1 pD(D2_7)
+1 pDis(B2_7)
+1 pG(G2_7)
-1 pLfrom(L3_23_7)
-1 pLfrom(L4_24_7)
-1 pLfrom(L5_25_7)
-1 pLto(L1_12_7)
= 0

c_e_KCL_real(2_8)_:
-1 pChar(B2_8)
-1 pD(D2_8)
+1 pDis(B2_8)
+1 pG(G2_8)
-1 pLfrom(L3_23_8)
-1 pLfrom(L4_24_8)
-1 pLfrom(L5_25_8)
-1 pLto(L1_12_8)
= 0

c_e_KCL_real(2_9)_:
-1 pChar(B2_9)
-1 pD(D2_9)
+1 pDis(B2_9)
+1 pG(G2_9)
-1 pLfrom(L3_23_9)
-1 pLfrom(L4_24_9)
-1 pLfrom(L5_25_9)
-1 pLto(L1_12_9)
= 0

c_e_KCL_real(3_0)_:
-1 pD(D3_0)
+1 pG(G3_0)
-1 pLfrom(L6_34_0)
-1 pLto(L3_23_0)
= 0

c_e_KCL_real(3_1)_:
-1 pD(D3_1)
+1 pG(G3_1)
-1 pLfrom(L6_34_1)
-1 pLto(L3_23_1)
= 0

c_e_KCL_real(3_2)_:
-1 pD(D3_2)
+1 pG(G3_2)
-1 pLfrom(L6_34_2)
-1 pLto(L3_23_2)
= 0

c_e_KCL_real(3_3)_:
-1 pD(D3_3)
+1 pG(G3_3)
-1 pLfrom(L6_34_3)
-1 pLto(L3_23_3)
= 0

c_e_KCL_real(3_4)_:
-1 pD(D3_4)
+1 pG(G3_4)
-1 pLfrom(L6_34_4)
-1 pLto(L3_23_4)
= 0

c_e_KCL_real(3_5)_:
-1 pD(D3_5)
+1 pG(G3_5)
-1 pLfrom(L6_34_5)
-1 pLto(L3_23_5)
= 0

c_e_KCL_real(3_6)_:
-1 pD(D3_6)
+1 pG(G3_6)
-1 pLfrom(L6_34_6)
-1 pLto(L3_23_6)
= 0

c_e_KCL_real(3_7)_:
-1 pD(D3_7)
+1 pG(G3_7)
-1 pLfrom(L6_34_7)
-1 pLto(L3_23_7)
= 0

c_e_KCL_real(3_8)_:
-1 pD(D3_8)
+1 pG(G3_8)
-1 pLfrom(L6_34_8)
-1 pLto(L3_23_8)
= 0

c_e_KCL_real(3_9)_:
-1 pD(D3_9)
+1 pG(G3_9)
-1 pLfrom(L6_34_9)
-1 pLto(L3_23_9)
= 0

c_e_KCL_real(4_0)_:
+1 pD(D4_0)
+1 pLfrom(L7_45_0)
+1 pLfromT(T1_47_0)
+1 pLfromT(T2_49_0)
+1 pLto(L4_24_0)
+1 pLto(L6_34_0)
= 0

c_e_KCL_real(4_1)_:
+1 pD(D4_1)
+1 pLfrom(L7_45_1)
+1 pLfromT(T1_47_1)
+1 pLfromT(T2_49_1)
+1 pLto(L4_24_1)
+1 pLto(L6_34_1)
= 0

c_e_KCL_real(4_2)_:
+1 pD(D4_2)
+1 pLfrom(L7_45_2)
+1 pLfromT(T1_47_2)
+1 pLfromT(T2_49_2)
+1 pLto(L4_24_2)
+1 pLto(L6_34_2)
= 0

c_e_KCL_real(4_3)_:
+1 pD(D4_3)
+1 pLfrom(L7_45_3)
+1 pLfromT(T1_47_3)
+1 pLfromT(T2_49_3)
+1 pLto(L4_24_3)
+1 pLto(L6_34_3)
= 0

c_e_KCL_real(4_4)_:
+1 pD(D4_4)
+1 pLfrom(L7_45_4)
+1 pLfromT(T1_47_4)
+1 pLfromT(T2_49_4)
+1 pLto(L4_24_4)
+1 pLto(L6_34_4)
= 0

c_e_KCL_real(4_5)_:
+1 pD(D4_5)
+1 pLfrom(L7_45_5)
+1 pLfromT(T1_47_5)
+1 pLfromT(T2_49_5)
+1 pLto(L4_24_5)
+1 pLto(L6_34_5)
= 0

c_e_KCL_real(4_6)_:
+1 pD(D4_6)
+1 pLfrom(L7_45_6)
+1 pLfromT(T1_47_6)
+1 pLfromT(T2_49_6)
+1 pLto(L4_24_6)
+1 pLto(L6_34_6)
= 0

c_e_KCL_real(4_7)_:
+1 pD(D4_7)
+1 pLfrom(L7_45_7)
+1 pLfromT(T1_47_7)
+1 pLfromT(T2_49_7)
+1 pLto(L4_24_7)
+1 pLto(L6_34_7)
= 0

c_e_KCL_real(4_8)_:
+1 pD(D4_8)
+1 pLfrom(L7_45_8)
+1 pLfromT(T1_47_8)
+1 pLfromT(T2_49_8)
+1 pLto(L4_24_8)
+1 pLto(L6_34_8)
= 0

c_e_KCL_real(4_9)_:
+1 pD(D4_9)
+1 pLfrom(L7_45_9)
+1 pLfromT(T1_47_9)
+1 pLfromT(T2_49_9)
+1 pLto(L4_24_9)
+1 pLto(L6_34_9)
= 0

c_e_KCL_real(5_0)_:
+1 pD(D5_0)
+1 pLfromT(T3_56_0)
+1 pLto(L2_15_0)
+1 pLto(L5_25_0)
+1 pLto(L7_45_0)
= 0

c_e_KCL_real(5_1)_:
+1 pD(D5_1)
+1 pLfromT(T3_56_1)
+1 pLto(L2_15_1)
+1 pLto(L5_25_1)
+1 pLto(L7_45_1)
= 0

c_e_KCL_real(5_2)_:
+1 pD(D5_2)
+1 pLfromT(T3_56_2)
+1 pLto(L2_15_2)
+1 pLto(L5_25_2)
+1 pLto(L7_45_2)
= 0

c_e_KCL_real(5_3)_:
+1 pD(D5_3)
+1 pLfromT(T3_56_3)
+1 pLto(L2_15_3)
+1 pLto(L5_25_3)
+1 pLto(L7_45_3)
= 0

c_e_KCL_real(5_4)_:
+1 pD(D5_4)
+1 pLfromT(T3_56_4)
+1 pLto(L2_15_4)
+1 pLto(L5_25_4)
+1 pLto(L7_45_4)
= 0

c_e_KCL_real(5_5)_:
+1 pD(D5_5)
+1 pLfromT(T3_56_5)
+1 pLto(L2_15_5)
+1 pLto(L5_25_5)
+1 pLto(L7_45_5)
= 0

c_e_KCL_real(5_6)_:
+1 pD(D5_6)
+1 pLfromT(T3_56_6)
+1 pLto(L2_15_6)
+1 pLto(L5_25_6)
+1 pLto(L7_45_6)
= 0

c_e_KCL_real(5_7)_:
+1 pD(D5_7)
+1 pLfromT(T3_56_7)
+1 pLto(L2_15_7)
+1 pLto(L5_25_7)
+1 pLto(L7_45_7)
= 0

c_e_KCL_real(5_8)_:
+1 pD(D5_8)
+1 pLfromT(T3_56_8)
+1 pLto(L2_15_8)
+1 pLto(L5_25_8)
+1 pLto(L7_45_8)
= 0

c_e_KCL_real(5_9)_:
+1 pD(D5_9)
+1 pLfromT(T3_56_9)
+1 pLto(L2_15_9)
+1 pLto(L5_25_9)
+1 pLto(L7_45_9)
= 0

c_e_KCL_real(6_0)_:
-1 pD(D6_0)
+1 pG(G4_0)
-1 pLfrom(L10_613_0)
-1 pLfrom(L8_611_0)
-1 pLfrom(L9_612_0)
-1 pLtoT(T3_56_0)
= 0

c_e_KCL_real(6_1)_:
-1 pD(D6_1)
+1 pG(G4_1)
-1 pLfrom(L10_613_1)
-1 pLfrom(L8_611_1)
-1 pLfrom(L9_612_1)
-1 pLtoT(T3_56_1)
= 0

c_e_KCL_real(6_2)_:
-1 pD(D6_2)
+1 pG(G4_2)
-1 pLfrom(L10_613_2)
-1 pLfrom(L8_611_2)
-1 pLfrom(L9_612_2)
-1 pLtoT(T3_56_2)
= 0

c_e_KCL_real(6_3)_:
-1 pD(D6_3)
+1 pG(G4_3)
-1 pLfrom(L10_613_3)
-1 pLfrom(L8_611_3)
-1 pLfrom(L9_612_3)
-1 pLtoT(T3_56_3)
= 0

c_e_KCL_real(6_4)_:
-1 pD(D6_4)
+1 pG(G4_4)
-1 pLfrom(L10_613_4)
-1 pLfrom(L8_611_4)
-1 pLfrom(L9_612_4)
-1 pLtoT(T3_56_4)
= 0

c_e_KCL_real(6_5)_:
-1 pD(D6_5)
+1 pG(G4_5)
-1 pLfrom(L10_613_5)
-1 pLfrom(L8_611_5)
-1 pLfrom(L9_612_5)
-1 pLtoT(T3_56_5)
= 0

c_e_KCL_real(6_6)_:
-1 pD(D6_6)
+1 pG(G4_6)
-1 pLfrom(L10_613_6)
-1 pLfrom(L8_611_6)
-1 pLfrom(L9_612_6)
-1 pLtoT(T3_56_6)
= 0

c_e_KCL_real(6_7)_:
-1 pD(D6_7)
+1 pG(G4_7)
-1 pLfrom(L10_613_7)
-1 pLfrom(L8_611_7)
-1 pLfrom(L9_612_7)
-1 pLtoT(T3_56_7)
= 0

c_e_KCL_real(6_8)_:
-1 pD(D6_8)
+1 pG(G4_8)
-1 pLfrom(L10_613_8)
-1 pLfrom(L8_611_8)
-1 pLfrom(L9_612_8)
-1 pLtoT(T3_56_8)
= 0

c_e_KCL_real(6_9)_:
-1 pD(D6_9)
+1 pG(G4_9)
-1 pLfrom(L10_613_9)
-1 pLfrom(L8_611_9)
-1 pLfrom(L9_612_9)
-1 pLtoT(T3_56_9)
= 0

c_e_KCL_real(7_0)_:
+1 pD(D7_0)
+1 pLfrom(L11_78_0)
+1 pLfrom(L12_79_0)
+1 pLtoT(T1_47_0)
= 0

c_e_KCL_real(7_1)_:
+1 pD(D7_1)
+1 pLfrom(L11_78_1)
+1 pLfrom(L12_79_1)
+1 pLtoT(T1_47_1)
= 0

c_e_KCL_real(7_2)_:
+1 pD(D7_2)
+1 pLfrom(L11_78_2)
+1 pLfrom(L12_79_2)
+1 pLtoT(T1_47_2)
= 0

c_e_KCL_real(7_3)_:
+1 pD(D7_3)
+1 pLfrom(L11_78_3)
+1 pLfrom(L12_79_3)
+1 pLtoT(T1_47_3)
= 0

c_e_KCL_real(7_4)_:
+1 pD(D7_4)
+1 pLfrom(L11_78_4)
+1 pLfrom(L12_79_4)
+1 pLtoT(T1_47_4)
= 0

c_e_KCL_real(7_5)_:
+1 pD(D7_5)
+1 pLfrom(L11_78_5)
+1 pLfrom(L12_79_5)
+1 pLtoT(T1_47_5)
= 0

c_e_KCL_real(7_6)_:
+1 pD(D7_6)
+1 pLfrom(L11_78_6)
+1 pLfrom(L12_79_6)
+1 pLtoT(T1_47_6)
= 0

c_e_KCL_real(7_7)_:
+1 pD(D7_7)
+1 pLfrom(L11_78_7)
+1 pLfrom(L12_79_7)
+1 pLtoT(T1_47_7)
= 0

c_e_KCL_real(7_8)_:
+1 pD(D7_8)
+1 pLfrom(L11_78_8)
+1 pLfrom(L12_79_8)
+1 pLtoT(T1_47_8)
= 0

c_e_KCL_real(7_9)_:
+1 pD(D7_9)
+1 pLfrom(L11_78_9)
+1 pLfrom(L12_79_9)
+1 pLtoT(T1_47_9)
= 0

c_e_KCL_real(8_0)_:
-1 pD(D8_0)
+1 pG(G5_0)
-1 pLto(L11_78_0)
= 0

c_e_KCL_real(8_1)_:
-1 pD(D8_1)
+1 pG(G5_1)
-1 pLto(L11_78_1)
= 0

c_e_KCL_real(8_2)_:
-1 pD(D8_2)
+1 pG(G5_2)
-1 pLto(L11_78_2)
= 0

c_e_KCL_real(8_3)_:
-1 pD(D8_3)
+1 pG(G5_3)
-1 pLto(L11_78_3)
= 0

c_e_KCL_real(8_4)_:
-1 pD(D8_4)
+1 pG(G5_4)
-1 pLto(L11_78_4)
= 0

c_e_KCL_real(8_5)_:
-1 pD(D8_5)
+1 pG(G5_5)
-1 pLto(L11_78_5)
= 0

c_e_KCL_real(8_6)_:
-1 pD(D8_6)
+1 pG(G5_6)
-1 pLto(L11_78_6)
= 0

c_e_KCL_real(8_7)_:
-1 pD(D8_7)
+1 pG(G5_7)
-1 pLto(L11_78_7)
= 0

c_e_KCL_real(8_8)_:
-1 pD(D8_8)
+1 pG(G5_8)
-1 pLto(L11_78_8)
= 0

c_e_KCL_real(8_9)_:
-1 pD(D8_9)
+1 pG(G5_9)
-1 pLto(L11_78_9)
= 0

c_e_KCL_real(9_0)_:
+1 pD(D9_0)
+1 pLfrom(L13_910_0)
+1 pLfrom(L14_914_0)
+1 pLto(L12_79_0)
+1 pLtoT(T2_49_0)
= 0

c_e_KCL_real(9_1)_:
+1 pD(D9_1)
+1 pLfrom(L13_910_1)
+1 pLfrom(L14_914_1)
+1 pLto(L12_79_1)
+1 pLtoT(T2_49_1)
= 0

c_e_KCL_real(9_2)_:
+1 pD(D9_2)
+1 pLfrom(L13_910_2)
+1 pLfrom(L14_914_2)
+1 pLto(L12_79_2)
+1 pLtoT(T2_49_2)
= 0

c_e_KCL_real(9_3)_:
+1 pD(D9_3)
+1 pLfrom(L13_910_3)
+1 pLfrom(L14_914_3)
+1 pLto(L12_79_3)
+1 pLtoT(T2_49_3)
= 0

c_e_KCL_real(9_4)_:
+1 pD(D9_4)
+1 pLfrom(L13_910_4)
+1 pLfrom(L14_914_4)
+1 pLto(L12_79_4)
+1 pLtoT(T2_49_4)
= 0

c_e_KCL_real(9_5)_:
+1 pD(D9_5)
+1 pLfrom(L13_910_5)
+1 pLfrom(L14_914_5)
+1 pLto(L12_79_5)
+1 pLtoT(T2_49_5)
= 0

c_e_KCL_real(9_6)_:
+1 pD(D9_6)
+1 pLfrom(L13_910_6)
+1 pLfrom(L14_914_6)
+1 pLto(L12_79_6)
+1 pLtoT(T2_49_6)
= 0

c_e_KCL_real(9_7)_:
+1 pD(D9_7)
+1 pLfrom(L13_910_7)
+1 pLfrom(L14_914_7)
+1 pLto(L12_79_7)
+1 pLtoT(T2_49_7)
= 0

c_e_KCL_real(9_8)_:
+1 pD(D9_8)
+1 pLfrom(L13_910_8)
+1 pLfrom(L14_914_8)
+1 pLto(L12_79_8)
+1 pLtoT(T2_49_8)
= 0

c_e_KCL_real(9_9)_:
+1 pD(D9_9)
+1 pLfrom(L13_910_9)
+1 pLfrom(L14_914_9)
+1 pLto(L12_79_9)
+1 pLtoT(T2_49_9)
= 0

c_e_KCL_real(10_0)_:
+1 pD(D10_0)
+1 pLfrom(L15_1011_0)
+1 pLto(L13_910_0)
= 0

c_e_KCL_real(10_1)_:
+1 pD(D10_1)
+1 pLfrom(L15_1011_1)
+1 pLto(L13_910_1)
= 0

c_e_KCL_real(10_2)_:
+1 pD(D10_2)
+1 pLfrom(L15_1011_2)
+1 pLto(L13_910_2)
= 0

c_e_KCL_real(10_3)_:
+1 pD(D10_3)
+1 pLfrom(L15_1011_3)
+1 pLto(L13_910_3)
= 0

c_e_KCL_real(10_4)_:
+1 pD(D10_4)
+1 pLfrom(L15_1011_4)
+1 pLto(L13_910_4)
= 0

c_e_KCL_real(10_5)_:
+1 pD(D10_5)
+1 pLfrom(L15_1011_5)
+1 pLto(L13_910_5)
= 0

c_e_KCL_real(10_6)_:
+1 pD(D10_6)
+1 pLfrom(L15_1011_6)
+1 pLto(L13_910_6)
= 0

c_e_KCL_real(10_7)_:
+1 pD(D10_7)
+1 pLfrom(L15_1011_7)
+1 pLto(L13_910_7)
= 0

c_e_KCL_real(10_8)_:
+1 pD(D10_8)
+1 pLfrom(L15_1011_8)
+1 pLto(L13_910_8)
= 0

c_e_KCL_real(10_9)_:
+1 pD(D10_9)
+1 pLfrom(L15_1011_9)
+1 pLto(L13_910_9)
= 0

c_e_KCL_real(11_0)_:
+1 pD(D11_0)
+1 pLto(L15_1011_0)
+1 pLto(L8_611_0)
= 0

c_e_KCL_real(11_1)_:
+1 pD(D11_1)
+1 pLto(L15_1011_1)
+1 pLto(L8_611_1)
= 0

c_e_KCL_real(11_2)_:
+1 pD(D11_2)
+1 pLto(L15_1011_2)
+1 pLto(L8_611_2)
= 0

c_e_KCL_real(11_3)_:
+1 pD(D11_3)
+1 pLto(L15_1011_3)
+1 pLto(L8_611_3)
= 0

c_e_KCL_real(11_4)_:
+1 pD(D11_4)
+1 pLto(L15_1011_4)
+1 pLto(L8_611_4)
= 0

c_e_KCL_real(11_5)_:
+1 pD(D11_5)
+1 pLto(L15_1011_5)
+1 pLto(L8_611_5)
= 0

c_e_KCL_real(11_6)_:
+1 pD(D11_6)
+1 pLto(L15_1011_6)
+1 pLto(L8_611_6)
= 0

c_e_KCL_real(11_7)_:
+1 pD(D11_7)
+1 pLto(L15_1011_7)
+1 pLto(L8_611_7)
= 0

c_e_KCL_real(11_8)_:
+1 pD(D11_8)
+1 pLto(L15_1011_8)
+1 pLto(L8_611_8)
= 0

c_e_KCL_real(11_9)_:
+1 pD(D11_9)
+1 pLto(L15_1011_9)
+1 pLto(L8_611_9)
= 0

c_e_KCL_real(12_0)_:
+1 pD(D12_0)
+1 pLfrom(L16_1213_0)
+1 pLto(L9_612_0)
= 0

c_e_KCL_real(12_1)_:
+1 pD(D12_1)
+1 pLfrom(L16_1213_1)
+1 pLto(L9_612_1)
= 0

c_e_KCL_real(12_2)_:
+1 pD(D12_2)
+1 pLfrom(L16_1213_2)
+1 pLto(L9_612_2)
= 0

c_e_KCL_real(12_3)_:
+1 pD(D12_3)
+1 pLfrom(L16_1213_3)
+1 pLto(L9_612_3)
= 0

c_e_KCL_real(12_4)_:
+1 pD(D12_4)
+1 pLfrom(L16_1213_4)
+1 pLto(L9_612_4)
= 0

c_e_KCL_real(12_5)_:
+1 pD(D12_5)
+1 pLfrom(L16_1213_5)
+1 pLto(L9_612_5)
= 0

c_e_KCL_real(12_6)_:
+1 pD(D12_6)
+1 pLfrom(L16_1213_6)
+1 pLto(L9_612_6)
= 0

c_e_KCL_real(12_7)_:
+1 pD(D12_7)
+1 pLfrom(L16_1213_7)
+1 pLto(L9_612_7)
= 0

c_e_KCL_real(12_8)_:
+1 pD(D12_8)
+1 pLfrom(L16_1213_8)
+1 pLto(L9_612_8)
= 0

c_e_KCL_real(12_9)_:
+1 pD(D12_9)
+1 pLfrom(L16_1213_9)
+1 pLto(L9_612_9)
= 0

c_e_KCL_real(13_0)_:
+1 pD(D13_0)
+1 pLfrom(L17_1314_0)
+1 pLto(L10_613_0)
+1 pLto(L16_1213_0)
= 0

c_e_KCL_real(13_1)_:
+1 pD(D13_1)
+1 pLfrom(L17_1314_1)
+1 pLto(L10_613_1)
+1 pLto(L16_1213_1)
= 0

c_e_KCL_real(13_2)_:
+1 pD(D13_2)
+1 pLfrom(L17_1314_2)
+1 pLto(L10_613_2)
+1 pLto(L16_1213_2)
= 0

c_e_KCL_real(13_3)_:
+1 pD(D13_3)
+1 pLfrom(L17_1314_3)
+1 pLto(L10_613_3)
+1 pLto(L16_1213_3)
= 0

c_e_KCL_real(13_4)_:
+1 pD(D13_4)
+1 pLfrom(L17_1314_4)
+1 pLto(L10_613_4)
+1 pLto(L16_1213_4)
= 0

c_e_KCL_real(13_5)_:
+1 pD(D13_5)
+1 pLfrom(L17_1314_5)
+1 pLto(L10_613_5)
+1 pLto(L16_1213_5)
= 0

c_e_KCL_real(13_6)_:
+1 pD(D13_6)
+1 pLfrom(L17_1314_6)
+1 pLto(L10_613_6)
+1 pLto(L16_1213_6)
= 0

c_e_KCL_real(13_7)_:
+1 pD(D13_7)
+1 pLfrom(L17_1314_7)
+1 pLto(L10_613_7)
+1 pLto(L16_1213_7)
= 0

c_e_KCL_real(13_8)_:
+1 pD(D13_8)
+1 pLfrom(L17_1314_8)
+1 pLto(L10_613_8)
+1 pLto(L16_1213_8)
= 0

c_e_KCL_real(13_9)_:
+1 pD(D13_9)
+1 pLfrom(L17_1314_9)
+1 pLto(L10_613_9)
+1 pLto(L16_1213_9)
= 0

c_e_KCL_real(14_0)_:
+1 pD(D14_0)
+1 pLto(L14_914_0)
+1 pLto(L17_1314_0)
= 0

c_e_KCL_real(14_1)_:
+1 pD(D14_1)
+1 pLto(L14_914_1)
+1 pLto(L17_1314_1)
= 0

c_e_KCL_real(14_2)_:
+1 pD(D14_2)
+1 pLto(L14_914_2)
+1 pLto(L17_1314_2)
= 0

c_e_KCL_real(14_3)_:
+1 pD(D14_3)
+1 pLto(L14_914_3)
+1 pLto(L17_1314_3)
= 0

c_e_KCL_real(14_4)_:
+1 pD(D14_4)
+1 pLto(L14_914_4)
+1 pLto(L17_1314_4)
= 0

c_e_KCL_real(14_5)_:
+1 pD(D14_5)
+1 pLto(L14_914_5)
+1 pLto(L17_1314_5)
= 0

c_e_KCL_real(14_6)_:
+1 pD(D14_6)
+1 pLto(L14_914_6)
+1 pLto(L17_1314_6)
= 0

c_e_KCL_real(14_7)_:
+1 pD(D14_7)
+1 pLto(L14_914_7)
+1 pLto(L17_1314_7)
= 0

c_e_KCL_real(14_8)_:
+1 pD(D14_8)
+1 pLto(L14_914_8)
+1 pLto(L17_1314_8)
= 0

c_e_KCL_real(14_9)_:
+1 pD(D14_9)
+1 pLto(L14_914_9)
+1 pLto(L17_1314_9)
= 0

c_e_KCL_reactive(1_0)_:
-1 qD(D1_0)
+1 qG(G1_0)
-1 qLfrom(L1_12_0)
-1 qLfrom(L2_15_0)
= 0

c_e_KCL_reactive(1_1)_:
-1 qD(D1_1)
+1 qG(G1_1)
-1 qLfrom(L1_12_1)
-1 qLfrom(L2_15_1)
= 0

c_e_KCL_reactive(1_2)_:
-1 qD(D1_2)
+1 qG(G1_2)
-1 qLfrom(L1_12_2)
-1 qLfrom(L2_15_2)
= 0

c_e_KCL_reactive(1_3)_:
-1 qD(D1_3)
+1 qG(G1_3)
-1 qLfrom(L1_12_3)
-1 qLfrom(L2_15_3)
= 0

c_e_KCL_reactive(1_4)_:
-1 qD(D1_4)
+1 qG(G1_4)
-1 qLfrom(L1_12_4)
-1 qLfrom(L2_15_4)
= 0

c_e_KCL_reactive(1_5)_:
-1 qD(D1_5)
+1 qG(G1_5)
-1 qLfrom(L1_12_5)
-1 qLfrom(L2_15_5)
= 0

c_e_KCL_reactive(1_6)_:
-1 qD(D1_6)
+1 qG(G1_6)
-1 qLfrom(L1_12_6)
-1 qLfrom(L2_15_6)
= 0

c_e_KCL_reactive(1_7)_:
-1 qD(D1_7)
+1 qG(G1_7)
-1 qLfrom(L1_12_7)
-1 qLfrom(L2_15_7)
= 0

c_e_KCL_reactive(1_8)_:
-1 qD(D1_8)
+1 qG(G1_8)
-1 qLfrom(L1_12_8)
-1 qLfrom(L2_15_8)
= 0

c_e_KCL_reactive(1_9)_:
-1 qD(D1_9)
+1 qG(G1_9)
-1 qLfrom(L1_12_9)
-1 qLfrom(L2_15_9)
= 0

c_e_KCL_reactive(2_0)_:
-1 qD(D2_0)
+1 qG(G2_0)
-1 qLfrom(L3_23_0)
-1 qLfrom(L4_24_0)
-1 qLfrom(L5_25_0)
-1 qLto(L1_12_0)
= 0

c_e_KCL_reactive(2_1)_:
-1 qD(D2_1)
+1 qG(G2_1)
-1 qLfrom(L3_23_1)
-1 qLfrom(L4_24_1)
-1 qLfrom(L5_25_1)
-1 qLto(L1_12_1)
= 0

c_e_KCL_reactive(2_2)_:
-1 qD(D2_2)
+1 qG(G2_2)
-1 qLfrom(L3_23_2)
-1 qLfrom(L4_24_2)
-1 qLfrom(L5_25_2)
-1 qLto(L1_12_2)
= 0

c_e_KCL_reactive(2_3)_:
-1 qD(D2_3)
+1 qG(G2_3)
-1 qLfrom(L3_23_3)
-1 qLfrom(L4_24_3)
-1 qLfrom(L5_25_3)
-1 qLto(L1_12_3)
= 0

c_e_KCL_reactive(2_4)_:
-1 qD(D2_4)
+1 qG(G2_4)
-1 qLfrom(L3_23_4)
-1 qLfrom(L4_24_4)
-1 qLfrom(L5_25_4)
-1 qLto(L1_12_4)
= 0

c_e_KCL_reactive(2_5)_:
-1 qD(D2_5)
+1 qG(G2_5)
-1 qLfrom(L3_23_5)
-1 qLfrom(L4_24_5)
-1 qLfrom(L5_25_5)
-1 qLto(L1_12_5)
= 0

c_e_KCL_reactive(2_6)_:
-1 qD(D2_6)
+1 qG(G2_6)
-1 qLfrom(L3_23_6)
-1 qLfrom(L4_24_6)
-1 qLfrom(L5_25_6)
-1 qLto(L1_12_6)
= 0

c_e_KCL_reactive(2_7)_:
-1 qD(D2_7)
+1 qG(G2_7)
-1 qLfrom(L3_23_7)
-1 qLfrom(L4_24_7)
-1 qLfrom(L5_25_7)
-1 qLto(L1_12_7)
= 0

c_e_KCL_reactive(2_8)_:
-1 qD(D2_8)
+1 qG(G2_8)
-1 qLfrom(L3_23_8)
-1 qLfrom(L4_24_8)
-1 qLfrom(L5_25_8)
-1 qLto(L1_12_8)
= 0

c_e_KCL_reactive(2_9)_:
-1 qD(D2_9)
+1 qG(G2_9)
-1 qLfrom(L3_23_9)
-1 qLfrom(L4_24_9)
-1 qLfrom(L5_25_9)
-1 qLto(L1_12_9)
= 0

c_e_KCL_reactive(3_0)_:
-1 qD(D3_0)
+1 qG(G3_0)
-1 qLfrom(L6_34_0)
-1 qLto(L3_23_0)
= 0

c_e_KCL_reactive(3_1)_:
-1 qD(D3_1)
+1 qG(G3_1)
-1 qLfrom(L6_34_1)
-1 qLto(L3_23_1)
= 0

c_e_KCL_reactive(3_2)_:
-1 qD(D3_2)
+1 qG(G3_2)
-1 qLfrom(L6_34_2)
-1 qLto(L3_23_2)
= 0

c_e_KCL_reactive(3_3)_:
-1 qD(D3_3)
+1 qG(G3_3)
-1 qLfrom(L6_34_3)
-1 qLto(L3_23_3)
= 0

c_e_KCL_reactive(3_4)_:
-1 qD(D3_4)
+1 qG(G3_4)
-1 qLfrom(L6_34_4)
-1 qLto(L3_23_4)
= 0

c_e_KCL_reactive(3_5)_:
-1 qD(D3_5)
+1 qG(G3_5)
-1 qLfrom(L6_34_5)
-1 qLto(L3_23_5)
= 0

c_e_KCL_reactive(3_6)_:
-1 qD(D3_6)
+1 qG(G3_6)
-1 qLfrom(L6_34_6)
-1 qLto(L3_23_6)
= 0

c_e_KCL_reactive(3_7)_:
-1 qD(D3_7)
+1 qG(G3_7)
-1 qLfrom(L6_34_7)
-1 qLto(L3_23_7)
= 0

c_e_KCL_reactive(3_8)_:
-1 qD(D3_8)
+1 qG(G3_8)
-1 qLfrom(L6_34_8)
-1 qLto(L3_23_8)
= 0

c_e_KCL_reactive(3_9)_:
-1 qD(D3_9)
+1 qG(G3_9)
-1 qLfrom(L6_34_9)
-1 qLto(L3_23_9)
= 0

c_e_KCL_reactive(4_0)_:
+1 qD(D4_0)
+1 qLfrom(L7_45_0)
+1 qLfromT(T1_47_0)
+1 qLfromT(T2_49_0)
+1 qLto(L4_24_0)
+1 qLto(L6_34_0)
= 0

c_e_KCL_reactive(4_1)_:
+1 qD(D4_1)
+1 qLfrom(L7_45_1)
+1 qLfromT(T1_47_1)
+1 qLfromT(T2_49_1)
+1 qLto(L4_24_1)
+1 qLto(L6_34_1)
= 0

c_e_KCL_reactive(4_2)_:
+1 qD(D4_2)
+1 qLfrom(L7_45_2)
+1 qLfromT(T1_47_2)
+1 qLfromT(T2_49_2)
+1 qLto(L4_24_2)
+1 qLto(L6_34_2)
= 0

c_e_KCL_reactive(4_3)_:
+1 qD(D4_3)
+1 qLfrom(L7_45_3)
+1 qLfromT(T1_47_3)
+1 qLfromT(T2_49_3)
+1 qLto(L4_24_3)
+1 qLto(L6_34_3)
= 0

c_e_KCL_reactive(4_4)_:
+1 qD(D4_4)
+1 qLfrom(L7_45_4)
+1 qLfromT(T1_47_4)
+1 qLfromT(T2_49_4)
+1 qLto(L4_24_4)
+1 qLto(L6_34_4)
= 0

c_e_KCL_reactive(4_5)_:
+1 qD(D4_5)
+1 qLfrom(L7_45_5)
+1 qLfromT(T1_47_5)
+1 qLfromT(T2_49_5)
+1 qLto(L4_24_5)
+1 qLto(L6_34_5)
= 0

c_e_KCL_reactive(4_6)_:
+1 qD(D4_6)
+1 qLfrom(L7_45_6)
+1 qLfromT(T1_47_6)
+1 qLfromT(T2_49_6)
+1 qLto(L4_24_6)
+1 qLto(L6_34_6)
= 0

c_e_KCL_reactive(4_7)_:
+1 qD(D4_7)
+1 qLfrom(L7_45_7)
+1 qLfromT(T1_47_7)
+1 qLfromT(T2_49_7)
+1 qLto(L4_24_7)
+1 qLto(L6_34_7)
= 0

c_e_KCL_reactive(4_8)_:
+1 qD(D4_8)
+1 qLfrom(L7_45_8)
+1 qLfromT(T1_47_8)
+1 qLfromT(T2_49_8)
+1 qLto(L4_24_8)
+1 qLto(L6_34_8)
= 0

c_e_KCL_reactive(4_9)_:
+1 qD(D4_9)
+1 qLfrom(L7_45_9)
+1 qLfromT(T1_47_9)
+1 qLfromT(T2_49_9)
+1 qLto(L4_24_9)
+1 qLto(L6_34_9)
= 0

c_e_KCL_reactive(5_0)_:
+1 qD(D5_0)
+1 qLfromT(T3_56_0)
+1 qLto(L2_15_0)
+1 qLto(L5_25_0)
+1 qLto(L7_45_0)
= 0

c_e_KCL_reactive(5_1)_:
+1 qD(D5_1)
+1 qLfromT(T3_56_1)
+1 qLto(L2_15_1)
+1 qLto(L5_25_1)
+1 qLto(L7_45_1)
= 0

c_e_KCL_reactive(5_2)_:
+1 qD(D5_2)
+1 qLfromT(T3_56_2)
+1 qLto(L2_15_2)
+1 qLto(L5_25_2)
+1 qLto(L7_45_2)
= 0

c_e_KCL_reactive(5_3)_:
+1 qD(D5_3)
+1 qLfromT(T3_56_3)
+1 qLto(L2_15_3)
+1 qLto(L5_25_3)
+1 qLto(L7_45_3)
= 0

c_e_KCL_reactive(5_4)_:
+1 qD(D5_4)
+1 qLfromT(T3_56_4)
+1 qLto(L2_15_4)
+1 qLto(L5_25_4)
+1 qLto(L7_45_4)
= 0

c_e_KCL_reactive(5_5)_:
+1 qD(D5_5)
+1 qLfromT(T3_56_5)
+1 qLto(L2_15_5)
+1 qLto(L5_25_5)
+1 qLto(L7_45_5)
= 0

c_e_KCL_reactive(5_6)_:
+1 qD(D5_6)
+1 qLfromT(T3_56_6)
+1 qLto(L2_15_6)
+1 qLto(L5_25_6)
+1 qLto(L7_45_6)
= 0

c_e_KCL_reactive(5_7)_:
+1 qD(D5_7)
+1 qLfromT(T3_56_7)
+1 qLto(L2_15_7)
+1 qLto(L5_25_7)
+1 qLto(L7_45_7)
= 0

c_e_KCL_reactive(5_8)_:
+1 qD(D5_8)
+1 qLfromT(T3_56_8)
+1 qLto(L2_15_8)
+1 qLto(L5_25_8)
+1 qLto(L7_45_8)
= 0

c_e_KCL_reactive(5_9)_:
+1 qD(D5_9)
+1 qLfromT(T3_56_9)
+1 qLto(L2_15_9)
+1 qLto(L5_25_9)
+1 qLto(L7_45_9)
= 0

c_e_KCL_reactive(6_0)_:
-1 qD(D6_0)
+1 qG(G4_0)
-1 qLfrom(L10_613_0)
-1 qLfrom(L8_611_0)
-1 qLfrom(L9_612_0)
-1 qLtoT(T3_56_0)
= 0

c_e_KCL_reactive(6_1)_:
-1 qD(D6_1)
+1 qG(G4_1)
-1 qLfrom(L10_613_1)
-1 qLfrom(L8_611_1)
-1 qLfrom(L9_612_1)
-1 qLtoT(T3_56_1)
= 0

c_e_KCL_reactive(6_2)_:
-1 qD(D6_2)
+1 qG(G4_2)
-1 qLfrom(L10_613_2)
-1 qLfrom(L8_611_2)
-1 qLfrom(L9_612_2)
-1 qLtoT(T3_56_2)
= 0

c_e_KCL_reactive(6_3)_:
-1 qD(D6_3)
+1 qG(G4_3)
-1 qLfrom(L10_613_3)
-1 qLfrom(L8_611_3)
-1 qLfrom(L9_612_3)
-1 qLtoT(T3_56_3)
= 0

c_e_KCL_reactive(6_4)_:
-1 qD(D6_4)
+1 qG(G4_4)
-1 qLfrom(L10_613_4)
-1 qLfrom(L8_611_4)
-1 qLfrom(L9_612_4)
-1 qLtoT(T3_56_4)
= 0

c_e_KCL_reactive(6_5)_:
-1 qD(D6_5)
+1 qG(G4_5)
-1 qLfrom(L10_613_5)
-1 qLfrom(L8_611_5)
-1 qLfrom(L9_612_5)
-1 qLtoT(T3_56_5)
= 0

c_e_KCL_reactive(6_6)_:
-1 qD(D6_6)
+1 qG(G4_6)
-1 qLfrom(L10_613_6)
-1 qLfrom(L8_611_6)
-1 qLfrom(L9_612_6)
-1 qLtoT(T3_56_6)
= 0

c_e_KCL_reactive(6_7)_:
-1 qD(D6_7)
+1 qG(G4_7)
-1 qLfrom(L10_613_7)
-1 qLfrom(L8_611_7)
-1 qLfrom(L9_612_7)
-1 qLtoT(T3_56_7)
= 0

c_e_KCL_reactive(6_8)_:
-1 qD(D6_8)
+1 qG(G4_8)
-1 qLfrom(L10_613_8)
-1 qLfrom(L8_611_8)
-1 qLfrom(L9_612_8)
-1 qLtoT(T3_56_8)
= 0

c_e_KCL_reactive(6_9)_:
-1 qD(D6_9)
+1 qG(G4_9)
-1 qLfrom(L10_613_9)
-1 qLfrom(L8_611_9)
-1 qLfrom(L9_612_9)
-1 qLtoT(T3_56_9)
= 0

c_e_KCL_reactive(7_0)_:
+1 qD(D7_0)
+1 qLfrom(L11_78_0)
+1 qLfrom(L12_79_0)
+1 qLtoT(T1_47_0)
= 0

c_e_KCL_reactive(7_1)_:
+1 qD(D7_1)
+1 qLfrom(L11_78_1)
+1 qLfrom(L12_79_1)
+1 qLtoT(T1_47_1)
= 0

c_e_KCL_reactive(7_2)_:
+1 qD(D7_2)
+1 qLfrom(L11_78_2)
+1 qLfrom(L12_79_2)
+1 qLtoT(T1_47_2)
= 0

c_e_KCL_reactive(7_3)_:
+1 qD(D7_3)
+1 qLfrom(L11_78_3)
+1 qLfrom(L12_79_3)
+1 qLtoT(T1_47_3)
= 0

c_e_KCL_reactive(7_4)_:
+1 qD(D7_4)
+1 qLfrom(L11_78_4)
+1 qLfrom(L12_79_4)
+1 qLtoT(T1_47_4)
= 0

c_e_KCL_reactive(7_5)_:
+1 qD(D7_5)
+1 qLfrom(L11_78_5)
+1 qLfrom(L12_79_5)
+1 qLtoT(T1_47_5)
= 0

c_e_KCL_reactive(7_6)_:
+1 qD(D7_6)
+1 qLfrom(L11_78_6)
+1 qLfrom(L12_79_6)
+1 qLtoT(T1_47_6)
= 0

c_e_KCL_reactive(7_7)_:
+1 qD(D7_7)
+1 qLfrom(L11_78_7)
+1 qLfrom(L12_79_7)
+1 qLtoT(T1_47_7)
= 0

c_e_KCL_reactive(7_8)_:
+1 qD(D7_8)
+1 qLfrom(L11_78_8)
+1 qLfrom(L12_79_8)
+1 qLtoT(T1_47_8)
= 0

c_e_KCL_reactive(7_9)_:
+1 qD(D7_9)
+1 qLfrom(L11_78_9)
+1 qLfrom(L12_79_9)
+1 qLtoT(T1_47_9)
= 0

c_e_KCL_reactive(8_0)_:
-1 qD(D8_0)
+1 qG(G5_0)
-1 qLto(L11_78_0)
= 0

c_e_KCL_reactive(8_1)_:
-1 qD(D8_1)
+1 qG(G5_1)
-1 qLto(L11_78_1)
= 0

c_e_KCL_reactive(8_2)_:
-1 qD(D8_2)
+1 qG(G5_2)
-1 qLto(L11_78_2)
= 0

c_e_KCL_reactive(8_3)_:
-1 qD(D8_3)
+1 qG(G5_3)
-1 qLto(L11_78_3)
= 0

c_e_KCL_reactive(8_4)_:
-1 qD(D8_4)
+1 qG(G5_4)
-1 qLto(L11_78_4)
= 0

c_e_KCL_reactive(8_5)_:
-1 qD(D8_5)
+1 qG(G5_5)
-1 qLto(L11_78_5)
= 0

