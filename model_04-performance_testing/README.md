This looks to test performance of a hillslope using a 2year period in 4 hours, to see who gets the furthest or done the fastest.

* baseline run is the standard set of options, running on 8 cores

* 1_no_ddivqdT turns _off_ the "Jacobian terms: d div q / dT" term
* 1_yes_ddivsurfqdT turns _on_ the "Jacobian terms: d div surface q / dT" term 
* 1_no_ddivqEdp turns _off_ the "Jacobian terms: d div K grad T / dp" term
* 1_yes_ddivhqdp turns _on_ the "Jacobian terms: d div hq / dp,T" term,
  and _off_ the subsurface PK's "supress advective terms in
  preconditioner" which is subsumed by the previous one
* 1_yes_surface_advection Turn _on_ surface advection term

* 2_surface_rel_perm changes to the "from above" strategy that was figured out in the intermediate scale model (NOT DONE)

* 4_boomer_2cycles does 3 cycle applications
* 4_boomer_3sweeps does 3 sweeps of the smoother
* 4_boomer_65strongthreshold sets the stron ghtreshold to 0.65
* 4_boomer_useblockindices sets the "use block indices" parameter for blocking
* 4_ilu Just use ILU instead and make sure boomer rocks still...

* 5_gmres_3its Don't oversolve linear solve
* 5_gmres_10its Don't oversolve linear solve part 2
* 5_none no linear solve

Results after 2 hours (or 2 years, whichever comes first)
3_manning_10		Cycle = 2566,  Time [days] = 729,  dt [days] = 1
3_manning_10-2		Cycle = 2566,  Time [days] = 729,  dt [days] = 1

4_ilu			Cycle = 2930,  Time [days] = 530.8218855022358,  dt [days] = 0.08905724888210427
4_ilu-2			Cycle = 2931,  Time [days] = 530.9109427511179,  dt [days] = 0.08905724888210427

1_no_ddivqEdp		Cycle = 2713,  Time [days] = 522.458836768365,  dt [days] = 0.005357768503707592
1_no_ddivqEdp-2		Cycle = 2528,  Time [days] = 521.4158428196291,  dt [days] = 0.008151417928618482

4_boomer_3sweeps	Cycle = 2765,  Time [days] = 522.0273838626541,  dt [days] = 0.002686306074503119
4_boomer_2cycles	Cycle = 2638,  Time [days] = 521.5306275694153,  dt [days] = 0.03725148363750325
4_boomer_65strongthresh Cycle = 2441,  Time [days] = 520.8005468547971,  dt [days] = 0.004492608986501755
1_yes_surface_advection Cycle = 2379,  Time [days] = 520.506102236168,  dt [days] = 0.002110673924908042
1_no_ddivqdT		Cycle = 2480,  Time [days] = 520.5986374145382,  dt [days] = 0.002033409692812711
5_gmres_10its		Cycle = 2351,  Time [days] = 520.3787545458547,  dt [days] = 0.006437024973608829

0_baseline		Cycle = 2312,  Time [days] = 520.2111012502723,  dt [days] = 0.002349341680226779
baseline-2		Cycle = 2511,  Time [days] = 521,  dt [days] = 0.004137618857677336

5_gmres_3its		Cycle = 2266,  Time [days] = 520.0445069071573,  dt [days] = 0.005787765759741888

1_yes_ddivsurfq		Cycle = 1515,  Time [days] = 256.8984375,  dt [days] = 0.01171875

1_yes_ddivhqdp		Cycle = 950,  Time [days] = 157.2946626768023,  dt [days] = 0.001701571941860796
2_surface_rel_perm
3_manning_0.1		Cycle = 2357,  Time [days] = 157.4276054625853,  dt [days] = 0.001216908437932876
4_boomer_useblockindice Cycle = 900,  Time [days] = 156.437683847581,  dt [days] = 0.001973639459034462
5_none			Cycle = 7345,  Time [days] = 155.485460855076,  dt [days] = 0.0001018229254124755


Repeats:
4_ilu
0_baseline
1_no_ddivqEdp
3_manning_10


Note: discharge plot shows that a manning coefficient of 10 is significantly different:
- slower to peak (by 1 day)
- slower to rising limb (by 1 day)
- higher peak (38%)
- less total annual discharge (15%)


This motivates a check for actual dependency on manning's n...:
3_manning_* are runs for 5 hours on 2 years (to make sure they finish
both years) to look at the dependency of discharge on Manning's n.



. print_result.sh | sort -k2 -n

4_boomer_3sweeps 3187
4_ilu 3359
4_boomer_2cycles 3513
1_no_ddivqEdp 3637
1_no_ddivqdT 3952
4_boomer_65strongthreshold 4090
5_gmres_10its 4134
0_baseline 4196
5_gmres_3its 4220
1_yes_surface_advection 4269
1_yes_ddivsurfqdT 5575
1_yes_ddivhqdp 12409
5_none 13682
4_boomer_useblockindices 13683