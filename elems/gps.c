/* gps.c - GPT element bindings */

/* This file has been generated automatically by 'makeelems' */

#include <stdlib.h>
#include "utils.h"
#include "gps.h"

/* Declarations of initialisation functions */
void catdog_init( struct inputlist *p ) ;
void catleg_init( struct inputlist *p ) ;
void onAxis_init( struct inputlist *p ) ;
void timesetfile_init( struct inputlist *p ) ;
void ionizer_init( struct inputlist *p ) ;
void multirmax_init( struct inputlist *p ) ;
void cathode_init( struct inputlist *p ) ;
void csol_init( struct inputlist *p ) ;
void kasten_init( struct inputlist *p ) ;
void myquad_init( struct inputlist *p ) ;
void dogleg_init( struct inputlist *p ) ;
void uturn_init( struct inputlist *p ) ;
void map3D3Sym_B_init( struct inputlist *p ) ;
void sectormagnet_fix_init( struct inputlist *p ) ;
void boxsolenoid_init( struct inputlist *p ) ;
void GBzmin_init( struct inputlist *p ) ;
void setyoffset_init( struct inputlist *p ) ;
void setxoffset_init( struct inputlist *p ) ;
void addzdiv_init( struct inputlist *p ) ;
void addydiv_init( struct inputlist *p ) ;
void addxdiv_init( struct inputlist *p ) ;
void setGByemittance_init( struct inputlist *p ) ;
void setGBxemittance_init( struct inputlist *p ) ;
void setGdist_init( struct inputlist *p ) ;
void setGBthetadist_init( struct inputlist *p ) ;
void setGBphidist_init( struct inputlist *p ) ;
void setGBrxydist_init( struct inputlist *p ) ;
void setGBzdist_init( struct inputlist *p ) ;
void setGBydist_init( struct inputlist *p ) ;
void setGBxdist_init( struct inputlist *p ) ;
void settdist_init( struct inputlist *p ) ;
void setrxydist_init( struct inputlist *p ) ;
void setphidist_init( struct inputlist *p ) ;
void setzdist_init( struct inputlist *p ) ;
void setydist_init( struct inputlist *p ) ;
void setxydistbmp_init( struct inputlist *p ) ;
void setxdist_init( struct inputlist *p ) ;
void setreduce_init( struct inputlist *p ) ;
void settcopy_init( struct inputlist *p ) ;
void setcurvature_init( struct inputlist *p ) ;
void setcopy_init( struct inputlist *p ) ;
void setmirror_init( struct inputlist *p ) ;
void setmove_init( struct inputlist *p ) ;
void setshuffle_init( struct inputlist *p ) ;
void settransform_init( struct inputlist *p ) ;
void setscale_init( struct inputlist *p ) ;
void setellipse_init( struct inputlist *p ) ;
void setfile_init( struct inputlist *p ) ;
void setstartxyzgrid_init( struct inputlist *p ) ;
void setstartline_init( struct inputlist *p ) ;
void setstartpar_init( struct inputlist *p ) ;
void setrmacrodist_init( struct inputlist *p ) ;
void setparticles_init( struct inputlist *p ) ;
void scattertorus_init( struct inputlist *p ) ;
void scattersphere_init( struct inputlist *p ) ;
void scatterplate_init( struct inputlist *p ) ;
void scatterpipe_init( struct inputlist *p ) ;
void scatteriris_init( struct inputlist *p ) ;
void scattercone_init( struct inputlist *p ) ;
void scatterbitmap_init( struct inputlist *p ) ;
void forwardscatter_init( struct inputlist *p ) ;
void copperscatter_init( struct inputlist *p ) ;
void spacecharge3Dtree_init( struct inputlist *p ) ;
void spacecharge2Dline_init( struct inputlist *p ) ;
void setcharge2Dcircle_init( struct inputlist *p ) ;
void spacecharge2Dcircle_init( struct inputlist *p ) ;
void spacecharge3Dmesh_init( struct inputlist *p ) ;
void spacecharge3Dclassic_init( struct inputlist *p ) ;
void spacecharge3D_init( struct inputlist *p ) ;
void settotalcharge_init( struct inputlist *p ) ;
void zminmax_init( struct inputlist *p ) ;
void xymax_init( struct inputlist *p ) ;
void stdxyzmax_init( struct inputlist *p ) ;
void rmax_init( struct inputlist *p ) ;
void multislit_init( struct inputlist *p ) ;
void Gminmax_init( struct inputlist *p ) ;
void tout_init( struct inputlist *p ) ;
void tcontinue_init( struct inputlist *p ) ;
void snapshot_init( struct inputlist *p ) ;
void screen_init( struct inputlist *p ) ;
void pp_init( struct inputlist *p ) ;
void outputvalue_init( struct inputlist *p ) ;
void dtmaxt_init( struct inputlist *p ) ;
void dtmaxpars_init( struct inputlist *p ) ;
void accuracy_init( struct inputlist *p ) ;
void acceptance_init( struct inputlist *p ) ;
void startcathode_init( struct inputlist *p ) ;
void setcathode_init( struct inputlist *p ) ;
void startpgb_init( struct inputlist *p ) ;
void startpar_init( struct inputlist *p ) ;
void startgrid_init( struct inputlist *p ) ;
void startcyl_init( struct inputlist *p ) ;
void start_init( struct inputlist *p ) ;
void setrmacro_init( struct inputlist *p ) ;
void ezcell_init( struct inputlist *p ) ;
void bz_init( struct inputlist *p ) ;
void bend_init( struct inputlist *p ) ;
void viscosity_init( struct inputlist *p ) ;
void setrlimit_init( struct inputlist *p ) ;
void run_init( struct inputlist *p ) ;
void randomize_init( struct inputlist *p ) ;
void gravity_init( struct inputlist *p ) ;
void drift_init( struct inputlist *p ) ;
void collision_init( struct inputlist *p ) ;
void ccsflip_init( struct inputlist *p ) ;
void ccs_init( struct inputlist *p ) ;
void map3D_Hcomplex_init( struct inputlist *p ) ;
void map3D_Ecomplex_init( struct inputlist *p ) ;
void map3D_TM_init( struct inputlist *p ) ;
void map3D_remove_init( struct inputlist *p ) ;
void map3D_V_init( struct inputlist *p ) ;
void map3D_B_init( struct inputlist *p ) ;
void map3D_E_init( struct inputlist *p ) ;
void map25D_TM_init( struct inputlist *p ) ;
void map25D_B_init( struct inputlist *p ) ;
void map25D_E_init( struct inputlist *p ) ;
void map2D_V_init( struct inputlist *p ) ;
void map2Dr_E_init( struct inputlist *p ) ;
void map2D_Et_init( struct inputlist *p ) ;
void map2D_E_init( struct inputlist *p ) ;
void map2D_B_init( struct inputlist *p ) ;
void map1D_TM_init( struct inputlist *p ) ;
void map1D_E_init( struct inputlist *p ) ;
void map1D_B_init( struct inputlist *p ) ;
void gauss00mf_init( struct inputlist *p ) ;
void unduplan_init( struct inputlist *p ) ;
void undueqfo_init( struct inputlist *p ) ;
void pointchargeset_init( struct inputlist *p ) ;
void pointcharge_init( struct inputlist *p ) ;
void platecharge_init( struct inputlist *p ) ;
void linecharge_init( struct inputlist *p ) ;
void erect_init( struct inputlist *p ) ;
void ehole_init( struct inputlist *p ) ;
void ecyl_init( struct inputlist *p ) ;
void circlecharge_init( struct inputlist *p ) ;
void solenoid_init( struct inputlist *p ) ;
void sextupole_init( struct inputlist *p ) ;
void sectormagnet_init( struct inputlist *p ) ;
void rectmagnet_init( struct inputlist *p ) ;
void rectcoil_init( struct inputlist *p ) ;
void quadrupole_init( struct inputlist *p ) ;
void magpoint_init( struct inputlist *p ) ;
void magplate_init( struct inputlist *p ) ;
void magline_init( struct inputlist *p ) ;
void magdipole_init( struct inputlist *p ) ;
void linecurrent_init( struct inputlist *p ) ;
void bzsolenoid_init( struct inputlist *p ) ;
void barmagnet_init( struct inputlist *p ) ;
void trwlinbm_init( struct inputlist *p ) ;
void trwlinac_init( struct inputlist *p ) ;
void trwcell_init( struct inputlist *p ) ;
void TMrectcavity_init( struct inputlist *p ) ;
void TM110cylcavity_init( struct inputlist *p ) ;
void TM010cylcavity_init( struct inputlist *p ) ;
void TErectcavity_init( struct inputlist *p ) ;

/* Keyword binding */
struct elementnames elems[] = {
  { "catdog" , catdog_init },
  { "catleg" , catleg_init },
  { "onAxis" , onAxis_init },
  { "timesetfile" , timesetfile_init },
  { "ionizer" , ionizer_init },
  { "multirmax" , multirmax_init },
  { "cathode" , cathode_init },
  { "csol" , csol_init },
  { "kasten" , kasten_init },
  { "myquad" , myquad_init },
  { "dogleg" , dogleg_init },
  { "uturn" , uturn_init },
  { "map3D3Sym_B" , map3D3Sym_B_init },
  { "sectormagnet_fix" , sectormagnet_fix_init },
  { "boxsolenoid" , boxsolenoid_init },
  { "GBzmin" , GBzmin_init },
  { "setyoffset" , setyoffset_init },
  { "setxoffset" , setxoffset_init },
  { "addzdiv" , addzdiv_init },
  { "addydiv" , addydiv_init },
  { "addxdiv" , addxdiv_init },
  { "setGByemittance" , setGByemittance_init },
  { "setGBxemittance" , setGBxemittance_init },
  { "setGdist" , setGdist_init },
  { "setGBthetadist" , setGBthetadist_init },
  { "setGBphidist" , setGBphidist_init },
  { "setGBrxydist" , setGBrxydist_init },
  { "setGBzdist" , setGBzdist_init },
  { "setGBydist" , setGBydist_init },
  { "setGBxdist" , setGBxdist_init },
  { "settdist" , settdist_init },
  { "setrxydist" , setrxydist_init },
  { "setphidist" , setphidist_init },
  { "setzdist" , setzdist_init },
  { "setydist" , setydist_init },
  { "setxydistbmp" , setxydistbmp_init },
  { "setxdist" , setxdist_init },
  { "setreduce" , setreduce_init },
  { "settcopy" , settcopy_init },
  { "setcurvature" , setcurvature_init },
  { "setcopy" , setcopy_init },
  { "setmirror" , setmirror_init },
  { "setmove" , setmove_init },
  { "setshuffle" , setshuffle_init },
  { "settransform" , settransform_init },
  { "setscale" , setscale_init },
  { "setellipse" , setellipse_init },
  { "setfile" , setfile_init },
  { "setstartxyzgrid" , setstartxyzgrid_init },
  { "setstartline" , setstartline_init },
  { "setstartpar" , setstartpar_init },
  { "setrmacrodist" , setrmacrodist_init },
  { "setparticles" , setparticles_init },
  { "scattertorus" , scattertorus_init },
  { "scattersphere" , scattersphere_init },
  { "scatterplate" , scatterplate_init },
  { "scatterpipe" , scatterpipe_init },
  { "scatteriris" , scatteriris_init },
  { "scattercone" , scattercone_init },
  { "scatterbitmap" , scatterbitmap_init },
  { "forwardscatter" , forwardscatter_init },
  { "copperscatter" , copperscatter_init },
  { "spacecharge3Dtree" , spacecharge3Dtree_init },
  { "spacecharge2Dline" , spacecharge2Dline_init },
  { "setcharge2Dcircle" , setcharge2Dcircle_init },
  { "spacecharge2Dcircle" , spacecharge2Dcircle_init },
  { "spacecharge3Dmesh" , spacecharge3Dmesh_init },
  { "spacecharge3Dclassic" , spacecharge3Dclassic_init },
  { "spacecharge3D" , spacecharge3D_init },
  { "settotalcharge" , settotalcharge_init },
  { "zminmax" , zminmax_init },
  { "xymax" , xymax_init },
  { "stdxyzmax" , stdxyzmax_init },
  { "rmax" , rmax_init },
  { "multislit" , multislit_init },
  { "Gminmax" , Gminmax_init },
  { "tout" , tout_init },
  { "tcontinue" , tcontinue_init },
  { "snapshot" , snapshot_init },
  { "screen" , screen_init },
  { "pp" , pp_init },
  { "outputvalue" , outputvalue_init },
  { "dtmaxt" , dtmaxt_init },
  { "dtmaxpars" , dtmaxpars_init },
  { "accuracy" , accuracy_init },
  { "acceptance" , acceptance_init },
  { "startcathode" , startcathode_init },
  { "setcathode" , setcathode_init },
  { "startpgb" , startpgb_init },
  { "startpar" , startpar_init },
  { "startgrid" , startgrid_init },
  { "startcyl" , startcyl_init },
  { "start" , start_init },
  { "setrmacro" , setrmacro_init },
  { "ezcell" , ezcell_init },
  { "bz" , bz_init },
  { "bend" , bend_init },
  { "viscosity" , viscosity_init },
  { "setrlimit" , setrlimit_init },
  { "run" , run_init },
  { "randomize" , randomize_init },
  { "gravity" , gravity_init },
  { "drift" , drift_init },
  { "collision" , collision_init },
  { "ccsflip" , ccsflip_init },
  { "ccs" , ccs_init },
  { "map3D_Hcomplex" , map3D_Hcomplex_init },
  { "map3D_Ecomplex" , map3D_Ecomplex_init },
  { "map3D_TM" , map3D_TM_init },
  { "map3D_remove" , map3D_remove_init },
  { "map3D_V" , map3D_V_init },
  { "map3D_B" , map3D_B_init },
  { "map3D_E" , map3D_E_init },
  { "map25D_TM" , map25D_TM_init },
  { "map25D_B" , map25D_B_init },
  { "map25D_E" , map25D_E_init },
  { "map2D_V" , map2D_V_init },
  { "map2Dr_E" , map2Dr_E_init },
  { "map2D_Et" , map2D_Et_init },
  { "map2D_E" , map2D_E_init },
  { "map2D_B" , map2D_B_init },
  { "map1D_TM" , map1D_TM_init },
  { "map1D_E" , map1D_E_init },
  { "map1D_B" , map1D_B_init },
  { "gauss00mf" , gauss00mf_init },
  { "unduplan" , unduplan_init },
  { "undueqfo" , undueqfo_init },
  { "pointchargeset" , pointchargeset_init },
  { "pointcharge" , pointcharge_init },
  { "platecharge" , platecharge_init },
  { "linecharge" , linecharge_init },
  { "erect" , erect_init },
  { "ehole" , ehole_init },
  { "ecyl" , ecyl_init },
  { "circlecharge" , circlecharge_init },
  { "solenoid" , solenoid_init },
  { "sextupole" , sextupole_init },
  { "sectormagnet" , sectormagnet_init },
  { "rectmagnet" , rectmagnet_init },
  { "rectcoil" , rectcoil_init },
  { "quadrupole" , quadrupole_init },
  { "magpoint" , magpoint_init },
  { "magplate" , magplate_init },
  { "magline" , magline_init },
  { "magdipole" , magdipole_init },
  { "linecurrent" , linecurrent_init },
  { "bzsolenoid" , bzsolenoid_init },
  { "barmagnet" , barmagnet_init },
  { "trwlinbm" , trwlinbm_init },
  { "trwlinac" , trwlinac_init },
  { "trwcell" , trwcell_init },
  { "TMrectcavity" , TMrectcavity_init },
  { "TM110cylcavity" , TM110cylcavity_init },
  { "TM010cylcavity" , TM010cylcavity_init },
  { "TErectcavity" , TErectcavity_init },
  { "", NULL },
} ;
