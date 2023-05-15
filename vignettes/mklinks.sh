# CHNOSZ/vignettes/mklinks.sh
# Add documentation links to vignettes
# 20190125 jmd

# anintro.html: add links to help topics
# Set background-image:none to remove underlines (from bootstrap theme)
sed -i 's/<code>?`CHNOSZ-package`<\/code>/<code><a href="..\/html\/CHNOSZ-package.html" style="background-image:none;">?`CHNOSZ-package`<\/a><\/code>/g' anintro.html
sed -i 's/<code>?basis<\/code>/<code><a href="..\/html\/basis.html" style="background-image:none;">?basis<\/a><\/code>/g' anintro.html
sed -i 's/<code>?mosaic<\/code>/<code><a href="..\/html\/mosaic.html" style="background-image:none;">?mosaic<\/a><\/code>/g' anintro.html
sed -i 's/<code>?buffer<\/code>/<code><a href="..\/html\/buffer.html" style="background-image:none;">?buffer<\/a><\/code>/g' anintro.html
sed -i 's/<code>?solubility<\/code>/<code><a href="..\/html\/solubility.html" style="background-image:none;">?solubility<\/a><\/code>/g' anintro.html
sed -i 's/<code>?ionize.aa<\/code>/<code><a href="..\/html\/ionize.aa.html" style="background-image:none;">?ionize.aa<\/a><\/code>/g' anintro.html
sed -i 's/<code>?count.aa<\/code>/<code><a href="..\/html\/util.fasta.html" style="background-image:none;">?count.aa<\/a><\/code>/g' anintro.html
sed -i 's/<code>?thermo<\/code>/<code><a href="..\/html\/thermo.html" style="background-image:none;">?thermo<\/a><\/code>/g' anintro.html
sed -i 's/<code>?hkf<\/code>/<code><a href="..\/html\/eos.html" style="background-image:none;">?hkf<\/a><\/code>/g' anintro.html
sed -i 's/<code>?cgl<\/code>/<code><a href="..\/html\/eos.html" style="background-image:none;">?cgl<\/a><\/code>/g' anintro.html
sed -i 's/<code>?water<\/code>/<code><a href="..\/html\/water.html" style="background-image:none;">?water<\/a><\/code>/g' anintro.html
sed -i 's/<code>?subcrt<\/code>/<code><a href="..\/html\/subcrt.html" style="background-image:none;">?subcrt<\/a><\/code>/g' anintro.html
sed -i 's/<code>?EOSregress<\/code>/<code><a href="..\/html\/EOSregress.html" style="background-image:none;">?EOSregress<\/a><\/code>/g' anintro.html
sed -i 's/<code>?taxonomy<\/code>/<code><a href="..\/html\/taxonomy.html" style="background-image:none;">?taxonomy<\/a><\/code>/g' anintro.html

# anintro.html: add links to function names
# Start at line 120 (below the TOC)
sed -i '312,$s/<code>info()<\/code>/<code><a href="..\/html\/info.html" style="background-image:none;">info()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>ZC()<\/code>/<code><a href="..\/html\/util.formula.html" style="background-image:none;">ZC()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>affinity()<\/code>/<code><a href="..\/html\/affinity.html" style="background-image:none;">affinity()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>thermo.refs()<\/code>/<code><a href="..\/html\/util.data.html" style="background-image:none;">thermo.refs()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>makeup()<\/code>/<code><a href="..\/html\/makeup.html" style="background-image:none;">makeup()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>as.chemical.formula()<\/code>/<code><a href="..\/html\/util.formula.html" style="background-image:none;">as.chemical.formula()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>subcrt()<\/code>/<code><a href="..\/html\/subcrt.html" style="background-image:none;">subcrt()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>T.units()<\/code>/<code><a href="..\/html\/util.units.html" style="background-image:none;">T.units()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>E.units()<\/code>/<code><a href="..\/html\/util.units.html" style="background-image:none;">E.units()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>P.units()<\/code>/<code><a href="..\/html\/util.units.html" style="background-image:none;">P.units()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>equilibrate()<\/code>/<code><a href="..\/html\/equilibrate.html" style="background-image:none;">equilibrate()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>diagram()<\/code>/<code><a href="..\/html\/diagram.html" style="background-image:none;">diagram()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>basis()<\/code>/<code><a href="..\/html\/basis.html" style="background-image:none;">basis()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>species()<\/code>/<code><a href="..\/html\/species.html" style="background-image:none;">species()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>reset()<\/code>/<code><a href="..\/html\/thermo.html" style="background-image:none;">reset()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>describe.reaction()<\/code>/<code><a href="..\/html\/util.expression.html" style="background-image:none;">describe.reaction()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>swap.basis()<\/code>/<code><a href="..\/html\/swap.basis.html" style="background-image:none;">swap.basis()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>water.lines()<\/code>/<code><a href="..\/html\/util.plot.html" style="background-image:none;">water.lines()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>ratlab()<\/code>/<code><a href="..\/html\/util.expression.html" style="background-image:none;">ratlab()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>mosaic()<\/code>/<code><a href="..\/html\/mosaic.html" style="background-image:none;">mosaic()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>convert()<\/code>/<code><a href="..\/html\/util.units.html" style="background-image:none;">convert()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>mod.buffer()<\/code>/<code><a href="..\/html\/buffer.html" style="background-image:none;">mod.buffer()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>reset\$/<code><a href="..\/html\/thermo.html" style="background-image:none;">reset<\/a>\$/g' anintro.html
sed -i '312,$s/<code>solubility()<\/code>/<code><a href="..\/html\/solubility.html" style="background-image:none;">solubility()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>ZC.col()<\/code>/<code><a href="..\/html\/util.plot.html" style="background-image:none;">ZC.col()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>aminoacids(\&quot;\&quot;)<\/code>/<code><a href="..\/html\/util.seq.html" style="background-image:none;">aminoacids(\&quot;\&quot;)<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>expr.species()<\/code>/<code><a href="..\/html\/util.expression.html" style="background-image:none;">expr.species()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>nonideal()<\/code>/<code><a href="..\/html\/nonideal.html" style="background-image:none;">nonideal()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>pinfo()<\/code>/<code><a href="..\/html\/protein.info.html" style="background-image:none;">pinfo()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>protein.length()<\/code>/<code><a href="..\/html\/protein.info.html" style="background-image:none;">protein.length()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>protein.formula()<\/code>/<code><a href="..\/html\/protein.info.html" style="background-image:none;">protein.formula()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>ionize.aa()<\/code>/<code><a href="..\/html\/ionize.aa.html" style="background-image:none;">ionize.aa()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>read.fasta()<\/code>/<code><a href="..\/html\/util.fasta.html" style="background-image:none;">read.fasta()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>protein.basis()<\/code>/<code><a href="..\/html\/protein.info.html" style="background-image:none;">protein.basis()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>unitize()<\/code>/<code><a href="..\/html\/util.misc.html" style="background-image:none;">unitize()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>axis.label()<\/code>/<code><a href="..\/html\/util.expression.html" style="background-image:none;">axis.label()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>seq2aa()<\/code>/<code><a href="..\/html\/add.protein.html" style="background-image:none;">seq2aa()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>add.protein()<\/code>/<code><a href="..\/html\/add.protein.html" style="background-image:none;">add.protein()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>aminoacids()<\/code>/<code><a href="..\/html\/util.seq.html" style="background-image:none;">aminoacids()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>add.OBIGT()<\/code>/<code><a href="..\/html\/add.OBIGT.html" style="background-image:none;">add.OBIGT()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>mod.OBIGT()<\/code>/<code><a href="..\/html\/add.OBIGT.html" style="background-image:none;">mod.OBIGT()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>checkGHS()<\/code>/<code><a href="..\/html\/util.data.html" style="background-image:none;">checkGHS()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>checkEOS()<\/code>/<code><a href="..\/html\/util.data.html" style="background-image:none;">checkEOS()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>equil.reaction()<\/code>/<code><a href="..\/html\/equilibrate.html" style="background-image:none;">equil.reaction()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>EOSregress()<\/code>/<code><a href="..\/html\/EOSregress.html" style="background-image:none;">EOSregress()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>RH2OBIGT()<\/code>/<code><a href="..\/html\/util.data.html" style="background-image:none;">RH2OBIGT()<\/a><\/code>/g' anintro.html
sed -i '312,$s/<code>retrieve()<\/code>/<code><a href="..\/html\/retrieve.html" style="background-image:none;">retrieve()<\/a><\/code>/g' anintro.html

# anintro.html: functions from R base packages
sed -i '312,$s/<code>install.packages/<code><a href="..\/..\/utils\/html\/install.packages.html" style="background-image:none;">install.packages<\/a>/g' anintro.html
sed -i '312,$s/<code>library/<code><a href="..\/..\/base\/html\/library.html" style="background-image:none;">library<\/a>/g' anintro.html
sed -i '312,$s/help()/<a href="..\/..\/utils\/html\/help.html" style="background-image:none;">help()<\/a>/g' anintro.html
sed -i '312,$s/help.start()/<a href="..\/..\/utils\/html\/help.start.html" style="background-image:none;">help.start()<\/a>/g' anintro.html
sed -i '312,$s/write.table/<a href="..\/..\/utils\/html\/write.table.html" style="background-image:none;">write.table<\/a>/g' anintro.html
sed -i '312,$s/<code>plot()<\/code>/<code><a href="..\/..\/graphics\/html\/plot.html" style="background-image:none;">plot()<\/a><\/code>/g' anintro.html
sed -i '312,$s/lapply()/<a href="..\/..\/base\/html\/lapply.html" style="background-image:none;">lapply()<\/a>/g' anintro.html
sed -i '312,$s/do.call()/<a href="..\/..\/base\/html\/do.call.html" style="background-image:none;">do.call()<\/a>/g' anintro.html
sed -i '312,$s/rbind()/<a href="..\/..\/base\/html\/cbind.html" style="background-image:none;">rbind()<\/a>/g' anintro.html
sed -i '312,$s/matplot()/<a href="..\/..\/graphics\/html\/matplot.html" style="background-image:none;">matplot()<\/a>/g' anintro.html
sed -i '312,$s/unlist()/<a href="..\/..\/base\/html\/unlist.html" style="background-image:none;">unlist()<\/a>/g' anintro.html
sed -i '312,$s/heat.colors()/<a href="..\/..\/grDevices\/html\/palettes.html" style="background-image:none;">heat.colors()<\/a>/g' anintro.html
sed -i '312,$s/read.csv()/<a href="..\/..\/utils\/html\/read.table.html" style="background-image:none;">read.csv()<\/a>/g' anintro.html
sed -i '312,$s/seq()/<a href="..\/..\/base\/html\/seq.html" style="background-image:none;">seq()<\/a>/g' anintro.html
sed -i '312,$s/substitute()/<a href="..\/..\/base\/html\/substitute.html" style="background-image:none;">substitute()<\/a>/g' anintro.html
sed -i '312,$s/for()/<a href="..\/..\/base\/html\/Control.html" style="background-image:none;">for()<\/a>/g' anintro.html
sed -i '312,$s/browseURL()/<a href="..\/..\/utils\/html\/browseURL.html" style="background-image:none;">browseURL()<\/a>/g' anintro.html
sed -i '312,$s/Sys.Date()/<a href="..\/..\/base\/html\/Sys.time.html" style="background-image:none;">Sys.Date()<\/a>/g' anintro.html

# Add links to OBIGT.html
sed -i 's/reset()/<a href="..\/html\/thermo.html" style="color:\ red;">reset()<\/a>/g' OBIGT.html
# Note closing ) so we don't clash with add.OBIGT(
sed -i 's/OBIGT()/<a href="..\/html\/thermo.html" style="color:\ red;">OBIGT()<\/a>/g' OBIGT.html
sed -i 's/thermo()/<a href="..\/html\/thermo.html" style="color:\ green;">thermo()<\/a>/g' OBIGT.html
sed -i 's/add.OBIGT(/<a href="..\/html\/add.OBIGT.html" style="color:\ red;">add.OBIGT<\/a>(/g' OBIGT.html
sed -i 's/water(/<a href="..\/html\/water.html" style="color:\ green;">water<\/a>(/g' OBIGT.html
sed -i 's/demo(/<a href="..\/demo">demo<\/a>(/g' OBIGT.html

# Add links to equilibrium.html 20200710
sed -i 's/equilibrate()/<a href="..\/html\/equilibrate.html">equilibrate()<\/a>/g' equilibrium.html
sed -i 's/solubility()/<a href="..\/html\/solubility.html">solubility()<\/a>/g' equilibrium.html
sed -i 's/add.OBIGT()/<a href="..\/html\/add.OBIGT.html">add.OBIGT()<\/a>/g' equilibrium.html

# Add links to multi-metal.html 20200716
sed -i 's/affinity()/<a href="..\/html\/affinity.html">affinity()<\/a>/g' multi-metal.html
sed -i 's/mash()/<a href="..\/html\/mix.html">mash()<\/a>/g' multi-metal.html
sed -i 's/diagram()/<a href="..\/html\/diagram.html">diagram()<\/a>/g' multi-metal.html
sed -i 's/mosaic()/<a href="..\/html\/mosaic.html">mosaic()<\/a>/g' multi-metal.html
sed -i 's/equilibrate()/<a href="..\/html\/equilibrate.html">equilibrate()<\/a>/g' multi-metal.html
sed -i 's/rebalance()/<a href="..\/html\/mix.html">rebalance()<\/a>/g' multi-metal.html
sed -i 's/ratlab()/<a href="..\/html\/util.expression.html">ratlab()<\/a>/g' multi-metal.html
sed -i 's/mix()/<a href="..\/html\/mix.html">mix()<\/a>/g' multi-metal.html
sed -i 's/mod.OBIGT()/<a href="..\/html\/add.OBIGT.html">mod.OBIGT()<\/a>/g' multi-metal.html
sed -i 's/retrieve()/<a href="..\/html\/retrieve.html">retrieve()<\/a>/g' multi-metal.html
sed -i 's/NaCl()/<a href="..\/html\/NaCl.html">NaCl()<\/a>/g' multi-metal.html
sed -i 's/solubility()/<a href="..\/html\/solubility.html">solubility()<\/a>/g' multi-metal.html
sed -i 's/subcrt()/<a href="..\/html\/subcrt.html">subcrt()<\/a>/g' multi-metal.html
sed -i 's/convert()/<a href="..\/html\/util.units.html">convert()<\/a>/g' multi-metal.html

# Add links to custom_data.Rmd 20230302
# Start at line 375 (below the TOC)
sed -i '375,$s/reset()/<a href="..\/html\/thermo.html" style="color:\ red;">reset()<\/a>/g' custom_data.html
# Note closing ) so we don't clash with add.OBIGT(, mod.OBIGT(, or logB.to.OBIGT(
sed -i '375,$s/>OBIGT()/><a href="..\/html\/thermo.html" style="color:\ red;">OBIGT()<\/a>/g' custom_data.html
sed -i '375,$s/add.OBIGT(/<a href="..\/html\/add.OBIGT.html" style="color:\ red;">add.OBIGT<\/a>(/g' custom_data.html
sed -i '375,$s/mod.OBIGT(/<a href="..\/html\/add.OBIGT.html" style="color:\ red;">mod.OBIGT<\/a>(/g' custom_data.html
sed -i '375,$s/logB.to.OBIGT(/<a href="..\/html\/logB.to.OBIGT.html" style="color:\ red;">logB.to.OBIGT<\/a>(/g' custom_data.html
sed -i '375,$s/basis(/<a href="..\/html\/basis.html" style="color:\ red;">basis<\/a>(/g' custom_data.html
sed -i '375,$s/species(/<a href="..\/html\/species.html" style="color:\ red;">species<\/a>(/g' custom_data.html
sed -i '375,$s/E.units(/<a href="..\/html\/util.units.html" style="color:\ red;">E.units<\/a>(/g' custom_data.html
sed -i '375,$s/info(/<a href="..\/html\/info.html" style="color:\ green;">info<\/a>(/g' custom_data.html
sed -i '375,$s/subcrt(/<a href="..\/html\/subcrt.html" style="color:\ green;">subcrt<\/a>(/g' custom_data.html
sed -i '375,$s/affinity(/<a href="..\/html\/affinity.html" style="color:\ green;">affinity<\/a>(/g' custom_data.html
sed -i '375,$s/thermo.refs(/<a href="..\/html\/util.data.html" style="color:\ green;">thermo.refs<\/a>(/g' custom_data.html
sed -i '375,$s/thermo(/<a href="..\/html\/thermo.html" style="color:\ green;">thermo<\/a>(/g' custom_data.html
sed -i '375,$s/checkGHS(/<a href="..\/html\/util.data.html" style="color:\ green;">checkGHS<\/a>(/g' custom_data.html
sed -i '375,$s/checkEOS(/<a href="..\/html\/util.data.html" style="color:\ green;">checkEOS<\/a>(/g' custom_data.html
