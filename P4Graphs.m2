needsPackage "TensorComplexes"
needsPackage "FourTiTwo"
--ring with variables for concentration and covariance matrix
--n = number vertices
--c = number colors
createRing = method()
createRing(ZZ) := (n) -> (
    c :=binomial(n+1,2);
    M := multiSubsets(toList(1..n),2);
    R := QQ[for k in M list s_(toSequence(k)),t_1..t_c,for k in M list p_(toSequence(k)) ];
    R
    );

sourceRing = method()
sourceRing(ZZ) := (n) -> (
    c :=binomial(n+1,2);
    M := multiSubsets(toList(1..n),2);
    R := QQ[for k in M list s_(toSequence(k))];
    R
    )

inverseIdeal = method()
inverseIdeal(ZZ,Ring,List) := (n,R,colors) -> (
    M := local M;
    k := local k;
    i := local i;
    l := local l;
    S := local S;
    shift := {1,1};
    S = QQ[for k from 0 to length(colors)-1 list t_k];
    K := map(S^n,S^n,{});
    for i from 0 to length(colors)-1 do(
	M = map(S^n, S^n, join(for k in colors_i list toSequence(k-shift)=>t_(i),for k in colors_i list toSequence(reverse(k-shift))=>t_(i)));
	K = K + M;
	);
    f := map(S,R,for x in  multiSubsets(toList(1..n),2) list (-1)^(x_0 + x_1)*determinant(submatrix'(K,{x_0-1}, {x_1-1})));
    kernel f
    )

--- Enter texQ I for a clean TeX output

texQ = method()
texQ Ideal := I -> (
    glist := flatten entries gens I;          -- list of generators
    m := #glist;
    clean := f -> (
        s := toString tex f;
        if (#s >= 2 and s#0 == "$" and s#(#s-1) == "$") then s = substring(s, 1, #s-2);
        s
    );
    if m == 0 then return "\\[\n\\begin{aligned}\n\\end{aligned}\n\\]\n";
    accstr := "\\[\n\\begin{aligned}\n";      -- accumulator (string)

    for k from 0 to m-1 do (
        sstr := clean (glist#k);              -- generator as a clean string
        if k < m-1 then
            accstr = concatenate{accstr, "&", sstr, ", \\\\\n"}
        else
            accstr = concatenate{accstr, "&", sstr, ".\n"}  -- final line gets a period
    );
    accstr = concatenate{accstr, "\\end{aligned}\n\\]\n"};
    accstr
)

------4-Path graph examples with 3 edges in the same color (PD)
------
n = 4 --  PD12: non-toric, Lie Algebra dimension 42 - got results after 2 hours
R =sourceRing(n) 
colors = {
  {{1,1},{2,2},{3,3}, {4,4}},
  {{1,2},{3,4},{2,3}}
}
Id12=inverseIdeal(n,R,colors)
netList(Id12_*)
dim Id12
texQ Id12
------
n = 4 --  PD11 non-toric, Lie Algebra dimension 8
R =sourceRing(n) 
colors = {  {{1,1}},
  {{2,2}},
  {{3,3}},
  {{4,4}},
  {{1,2},{3,4},{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PD10 toric!, Lie Algebra dimension 14
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2},{4,4}},
  {{3,3}},
  {{1,2},{3,4},{2,3}}
}
Id10=inverseIdeal(n,R,colors)
netList(Id10_*)
dim Id10
texQ Id10
------
n = 4 --  PD9 non-toric, Lie Algebra dimension 3
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2},{3,3}},
  {{4,4}},
  {{1,2},{3,4},{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PD8 non-toric, Lie Algebra dimension 14
R =sourceRing(n) 
colors = {
  {{1,1},{2,2}},
  {{3,3}},
  {{4,4}},
  {{1,2},{3,4},{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PD7: toric!, Lie Algebra dimension 43 took forever to check if ideal is bionomial after the coordinate change
R =sourceRing(n) 
colors = {
  {{1,1},{3,3}},
  {{2,2},{4,4}},
  {{1,2},{3,4},{2,3}}
}
Id7=inverseIdeal(n,R,colors)
netList(Id7_*)
dim Id7
texQ Id7
------
n = 4 --  PD6: non-toric, Lie Algebra dimension 43 (took 20+ hours to diagonalizing torus)
R =sourceRing(n) 
colors = {
  {{1,1},{4,4}},
  {{2,2},{3,3}},
  {{1,2},{3,4},{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PD5: toric!, Lie Algebra dimension 25
R =sourceRing(n) 
colors = {
  {{1,1},{2,2}},
  {{3,3},{4,4}},
  {{1,2},{3,4},{2,3}}
}
Id5=inverseIdeal(n,R,colors)
netList(Id5_*)
dim Id5
texQ Id5
------
n = 4 --  PD4 (the same as PD1): toric!, Lie Algebra dimension 37 (the same as PD1)
R =sourceRing(n) 
colors = {
  {{4,4}},
  {{1,1},{2,2}, {3,3}},
  {{1,2},{3,4},{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PD3 (the same as PD2): toric!, Lie Algebra dimension 37 
R =sourceRing(n) 
colors = {
  {{3,3}},
  {{1,1},{2,2}, {4,4}},
  {{1,2},{3,4},{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PD2: toric!, Lie Algebra dimension 37
R =sourceRing(n) 
colors = {
  {{2,2}},
  {{1,1},{3,3}, {4,4}},
  {{1,2},{3,4},{2,3}}
}
Id2=inverseIdeal(n,R,colors)
netList(Id2_*)
dim Id2
texQ Id2
------
n = 4 --  PD1: toric!, Lie Algebra dimension 37
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2},{3,3}, {4,4}},
  {{1,2},{3,4},{2,3}}
}
Id1=inverseIdeal(n,R,colors)
netList(Id1_*)
dim Id1
texQ Id1

------4-Path graph examples with 3 edges colored distinctly (PC)
------
n = 4 --  PC12: non-toric, Lie Algebra dimension 12
R =sourceRing(n) 
colors = {
  {{1,1},{2,2},{3,3}, {4,4}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PC11 toric!, Lie Algebra dimension 8
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2}},
  {{3,3}},
  {{4,4}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
Ic11=inverseIdeal(n,R,colors)
netList(Ic11_*)
dim Ic11
texQ Ic11
------
n = 4 --  PC10 non-toric, Lie Algebra dimension 4
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2},{4,4}},
  {{3,3}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PC9 non-toric, Lie Algebra dimension 3
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2},{3,3}},
  {{4,4}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PC8 toric!, Lie Algebra dimension 9=8
R =sourceRing(n) 
colors = {
  {{1,1},{2,2}},
  {{3,3}},
  {{4,4}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
Ic8=inverseIdeal(n,R,colors)
netList(Ic8_*)
dim Ic8
texQ Ic8
-------
n = 4 --  PC7: non-toric, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{1,1},{3,3}},
  {{2,2}, {4,4}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PC6: non-toric, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{1,1},{4,4}},
  {{2,2}, {3,3}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PC5: non-toric, Lie Algebra dimension 3
R =sourceRing(n) 
colors = {
  {{1,1},{2,2}},
  {{3,3}, {4,4}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PC4 (Same with PC1): non-toric, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{4,4}},
  {{1,1}, {2,2}, {3,3}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PC3 (Same with PC2): non-toric, Lie Algebra dimension 4 
R =sourceRing(n) 
colors = {
  {{3,3}},
  {{1,1}, {2,2}, {4,4}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PC2: non-toric, Lie Algebra dimension 4
R =sourceRing(n) 
colors = {
  {{2,2}},
  {{1,1}, {3,3}, {4,4}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PC1: non-toric, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2}, {3,3}, {4,4}},
  {{1,2}},
  {{3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------4-Path graph examples with 2 edge colors: 12&34; 12
------
n = 4 --  PB12: non-toric, Lie Algebra dimension 42 (computation time 1h+)
R =sourceRing(n) 
colors = {
  {{1,1},{2,2},{3,3}, {4,4}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB11 non-toric, Lie Algebra dimension 8
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2}},
  {{3,3}},
  {{4,4}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB10 non-toric, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2},{4,4}},
  {{3,3}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB9 non-toric, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2},{3,3}},
  {{4,4}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB8 non-toric, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{1,1},{2,2}},
  {{3,3}},
  {{4,4}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB7: non-toric, Lie Algebra dimension 13
R =sourceRing(n) 
colors = {
  {{1,1}, {3,3}},
  {{2,2}, {4,4}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB6: toric!, Lie Algebra dimension 53 (know it through RCOP, but it took forever to check if ideal is bionomial after coordinate change)
R =sourceRing(n) 
colors = {
  {{1,1},{4,4}},
  {{2,2}, {3,3}},
  {{1,2}, {3,4}},
  {{2,3}}
}
Ib6=inverseIdeal(n,R,colors)
netList(Ib6_*)
dim Ib6
texQ Ib6
------
n = 4 --  PB5: non-toric, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{1,1},{2,2}},
  {{3,3}, {4,4}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB4: non-toric, Lie Algebra dimension 14
R =sourceRing(n) 
colors = {
  {{4,4}},
  {{1,1}, {2,2}, {3,3}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB3: non-toric, Lie Algebra dimension 14
R =sourceRing(n) 
colors = {
  {{3,3}},
  {{1,1}, {2,2}, {4,4}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB2: non-toric, Lie Algebra dimension 14
R =sourceRing(n) 
colors = {
  {{2,2}},
  {{1,1}, {3,3}, {4,4}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PB1: non-toric, Lie Algebra dimension 1
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2}, {3,3}, {4,4}},
  {{1,2}, {3,4}},
  {{2,3}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I

------PA: P4 graphs with 2 edge colors: 12&23; 34

------
n = 4 --  PA12 non-toric exmaple 8, Lie Algebra dimension 25
R =sourceRing(n) 
colors = {
  {{1,1},{2,2},{3,3},{4,4}},
  {{1,2}, {2,3}},
  {{3,4}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PA11 toric!, Lie Algebra dimension 8
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2}},
  {{3,3}},
  {{4,4}},
  {{1,2}, {2,3}},
  {{3,4}}
}
Ia11=inverseIdeal(n,R,colors)
netList(Ia11_*)
dim Ia11
texQ Ia11
------
n = 4 --  PA10 non-toric, Lie Algebra dimension 3
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2},{4,4}},
  {{3,3}},
  {{1,2}, {2,3}},
  {{3,4}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PA9 non-toric, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{1,1}},
  {{2,2},{3,3}},
  {{4,4}},
  {{1,2}, {2,3}},
  {{3,4}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PA8 toric!, Lie Algebra dimension 21
R =sourceRing(n) 
colors = {
  {{1,1},{2,2}},
  {{3,3}},
  {{4,4}},
  {{1,2}, {2,3}},
  {{3,4}}
}
Ia8=inverseIdeal(n,R,colors)
netList(Ia8_*)
dim Ia8
texQ Ia8
------
n = 4 -- PA7 non-toric exmaple 6, Lie Algebra dimension 2
R =sourceRing(n) 
colors = {
  {{1,1},{3,3}},
  {{2,2},{4,4}}, 
  {{1,2}, {2,3}},
  {{3,4}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 -- PA6 non-toric exmaple 5, Lie Algebra dimension 14
R =sourceRing(n) 
colors = {
  {{2,2}, {3,3}}, 
  {{1,1},{4,4}},
  {{1,2}, {2,3}},
  {{3,4}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 -- PA5  toric!, Lie Algebra dimension 14
R =sourceRing(n) 
colors = {
  {{1,1},{2,2}}, 
  {{3,3},{4,4}},
  {{1,2}, {2,3}},
  {{3,4}}
}
Ia5=inverseIdeal(n,R,colors)
netList(Ia5_*)
dim Ia5
texQ Ia5
------
n = 4 -- PA4  non-toric exmaple 4, Lie Algebra dimension 14
R =sourceRing(n) 
colors = {
  {{1,1},{2,2}, {3,3}}, 
  {{4,4}},
  {{1,2}, {2,3}},
  {{3,4}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 --  PA3 non-toric exmaple 7, Lie Algebra dimension 15
R =sourceRing(n) 
colors = {
  {{1,1},{2,2},{4,4}},
  {{3,3}}, 
  {{1,2}, {2,3}},
  {{3,4}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
texQ I
------
n = 4 -- PA2 non-toric exmaple 1, Lie Algebra dimension 14
R =sourceRing(n) 
colors = {
  {{1,1}, {3,3}, {4,4}},
  {{2,2}},
  {{1,2}, {2,3}},
  {{3,4}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I
------
n = 4 -- PA1 non-toric exmaple 2, Lie Algebra dimension 3
R =sourceRing(n) 
colors = {
  {{2,2}, {3,3}, {4,4}},
  {{1,1}},
  {{1,2}, {2,3}},
  {{3,4}}
}
I=inverseIdeal(n,R,colors)
netList(I_*)
dim I

