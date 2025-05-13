# CHNOSZ/vignettes/postprocess.sh
# Add documentation links to vignettes
# 20190125 jmd

### Start processing anintro.html

# Highlight NOTE in all caps
sed -i 's/NOTE/<span class="highlight">NOTE<\/span>/g' anintro.html

# Functions without side effects (green)
# Set background-image:none to remove underlines (from bootstrap theme)
sed -i 's/<code>info()/<code><a href="..\/html\/info.html" style="background-image:none;color:green;">info()<\/a>/g' anintro.html
sed -i 's/<code>retrieve()/<code><a href="..\/html\/retrieve.html" style="background-image:none;color:green;">retrieve()<\/a>/g' anintro.html
sed -i 's/<code>subcrt()/<code><a href="..\/html\/subcrt.html" style="background-image:none;color:green;">subcrt()<\/a>/g' anintro.html
sed -i 's/<code>affinity()/<code><a href="..\/html\/affinity.html" style="background-image:none;color:green;">affinity()<\/a>/g' anintro.html
sed -i 's/<code>mosaic()/<code><a href="..\/html\/mosaic.html" style="background-image:none;color:green;">mosaic()<\/a>/g' anintro.html
sed -i 's/<code>equilibrate()/<code><a href="..\/html\/equilibrate.html" style="background-image:none;color:green;">equilibrate()<\/a>/g' anintro.html
sed -i 's/<code>solubility()/<code><a href="..\/html\/solubility.html" style="background-image:none;color:green;">solubility()<\/a>/g' anintro.html
sed -i 's/<code>diagram()/<code><a href="..\/html\/diagram.html" style="background-image:none;color:green;">diagram()<\/a>/g' anintro.html
sed -i 's/<code>NaCl()/<code><a href="..\/html\/NaCl.html" style="background-image:none;color:green;">NaCl()<\/a>/g' anintro.html
sed -i 's/<code>affinity(return.buffer\ =\ TRUE)/<code><a href="..\/html\/affinity.html" style="background-image:none;color:green;">affinity(return.buffer\ =\ TRUE)<\/a>/g' anintro.html
sed -i 's/<code>ionize.aa()/<code><a href="..\/html\/ionize.aa.html" style="background-image:none;color:green;">ionize.aa()<\/a>/g' anintro.html
sed -i 's/<code>equilibrate(loga.balance\ =\ 0)/<code><a href="..\/html\/equilibrate.html" style="background-image:none;color:green;">equilibrate(loga.balance\ =\ 0)<\/a>/g' anintro.html
sed -i 's/<code>makeup()/<code><a href="..\/html\/makeup.html" style="background-image:none;color:green;">makeup()<\/a>/g' anintro.html
sed -i 's/<code>Berman()/<code><a href="..\/html\/Berman.html" style="background-image:none;color:green;">Berman()<\/a>/g' anintro.html

# Functions with side effects (red)
sed -i 's/<code>basis()/<code><a href="..\/html\/basis.html" style="background-image:none;color:red;">basis()<\/a>/g' anintro.html
sed -i 's/<code>swap.basis()/<code><a href="..\/html\/swap.basis.html" style="background-image:none;color:red;">swap.basis()<\/a>/g' anintro.html
sed -i 's/<code>species()/<code><a href="..\/html\/species.html" style="background-image:none;color:red;">species()<\/a>/g' anintro.html
sed -i 's/<code>OBIGT()/<code><a href="..\/html\/OBIGT.html" style="background-image:none;color:red;">OBIGT()<\/a>/g' anintro.html
sed -i 's/<code>add.OBIGT()/<code><a href="..\/html\/add.OBIGT.html" style="background-image:none;color:red;">add.OBIGT()<\/a>/g' anintro.html
sed -i 's/<code>logK.to.OBIGT()/<code><a href="..\/html\/logK.to.OBIGT.html" style="background-image:none;color:red;">logK.to.OBIGT()<\/a>/g' anintro.html
sed -i 's/<code>reset()/<code><a href="..\/html\/reset.html" style="background-image:none;color:red;">reset()<\/a>/g' anintro.html
sed -i 's/<code>thermo()/<code><a href="..\/html\/thermo.html" style="background-image:none;color:red;">thermo()<\/a>/g' anintro.html
sed -i 's/<code>mod.buffer()/<code><a href="..\/html\/mod.buffer.html" style="background-image:none;color:red;">mod.buffer()<\/a>/g' anintro.html
sed -i 's/<code>add.protein()/<code><a href="..\/html\/add.protein.html" style="background-image:none;color:red;">add.protein()<\/a>/g' anintro.html
sed -i 's/<code>nonideal()/<code><a href="..\/html\/nonideal.html" style="background-image:none;color:red;">nonideal()<\/a>/g' anintro.html

# Functions with different target page from function name
sed -i 's/<code>demos()/<code><a href="..\/html\/examples.html" style="background-image:none;color:green;">demos()<\/a>/g' anintro.html
sed -i 's/<code>convert()/<code><a href="..\/html\/util.units.html" style="background-image:none;color:green;">convert()<\/a>/g' anintro.html
sed -i 's/<code>water.lines()/<code><a href="..\/html\/util.plot.html" style="background-image:none;color:green;">water.lines()<\/a>/g' anintro.html
sed -i 's/<code>axis.label()/<code><a href="..\/html\/util.expression.html" style="background-image:none;color:green;">axis.label()<\/a>/g' anintro.html
sed -i 's/<code>expr.species()/<code><a href="..\/html\/util.expression.html" style="background-image:none;color:green;">expr.species()<\/a>/g' anintro.html
sed -i 's/<code>describe.reaction()/<code><a href="..\/html\/util.expression.html" style="background-image:none;color:green;">describe.reaction()<\/a>/g' anintro.html
sed -i 's/<code>describe.property()/<code><a href="..\/html\/util.expression.html" style="background-image:none;color:green;">describe.property()<\/a>/g' anintro.html
sed -i 's/<code>describe.basis()/<code><a href="..\/html\/util.expression.html" style="background-image:none;color:green;">describe.basis()<\/a>/g' anintro.html
sed -i 's/<code>lT()/<code><a href="..\/html\/util.legend.html" style="background-image:none;color:green;">lT()<\/a>/g' anintro.html
sed -i 's/<code>lP()/<code><a href="..\/html\/util.legend.html" style="background-image:none;color:green;">lP()<\/a>/g' anintro.html
sed -i 's/<code>lTP()/<code><a href="..\/html\/util.legend.html" style="background-image:none;color:green;">lTP()<\/a>/g' anintro.html
sed -i 's/<code>thermo.refs()/<code><a href="..\/html\/util.data.html" style="background-image:none;color:green;">thermo.refs()<\/a>/g' anintro.html
sed -i 's/<code>mod.OBIGT()/<code><a href="..\/html\/add.OBIGT.html" style="background-image:none;color:red;">mod.OBIGT()<\/a>/g' anintro.html
sed -i 's/<code>T.units()/<code><a href="..\/html\/util.units.html" style="background-image:none;color:red;">T.units()<\/a>/g' anintro.html
sed -i 's/<code>P.units()/<code><a href="..\/html\/util.units.html" style="background-image:none;color:red;">P.units()<\/a>/g' anintro.html
sed -i 's/<code>E.units()/<code><a href="..\/html\/util.units.html" style="background-image:none;color:red;">E.units()<\/a>/g' anintro.html
sed -i 's/<code>ZC()/<code><a href="..\/html\/util.formula.html" style="background-image:none;color:green;">ZC()<\/a>/g' anintro.html
sed -i 's/<code>unitize()/<code><a href="..\/html\/util.misc.html" style="background-image:none;color:green;">unitize()<\/a>/g' anintro.html
sed -i 's/<code>as.chemical.formula()/<code><a href="..\/html\/util.formula.html" style="background-image:none;color:green;">as.chemical.formula()<\/a>/g' anintro.html

# Help topics
sed -i 's/?info/<a href="..\/html\/info.html" style="background-image:none;color:green;">?info<\/a>/g' anintro.html
sed -i 's/?subcrt/<a href="..\/html\/subcrt.html" style="background-image:none;color:green;">?subcrt<\/a>/g' anintro.html

# Functions from R distribution packages
sed -i 's/<code>library/<code><a href="..\/..\/base\/html\/library.html">library<\/a>/g' anintro.html
sed -i 's/<code>bquote()/<code><a href="..\/..\/base\/html\/bquote.html">bquote()<\/a>/g' anintro.html
sed -i 's/<code>plotmath()/<code><a href="..\/..\/grDevices\/html\/plotmath.html">plotmath()<\/a>/g' anintro.html
sed -i 's/<code>install.packages/<code><a href="..\/..\/utils\/html\/install.packages.html">install.packages<\/a>/g' anintro.html
sed -i 's/<code>help(/<code><a href="..\/..\/utils\/html\/help.html">help(<\/a>/g' anintro.html
sed -i 's/<code>help.start/<code><a href="..\/..\/utils\/html\/help.start.html">help.start<\/a>/g' anintro.html
sed -i 's/<code>maintainer/<code><a href="..\/..\/utils\/html\/maintainer.html">maintainer<\/a>/g' anintro.html
sed -i 's/<code>as.expression/<code><a href="..\/..\/base\/html\/expression.html">as.expression<\/a>/g' anintro.html
sed -i 's/<code>legend/<code><a href="..\/..\/graphics\/html\/legend.html">legend<\/a>/g' anintro.html
sed -i 's/<code>log10()/<code><a href="..\/..\/base\/html\/Log.html">log10()<\/a>/g' anintro.html

### End processing anintro.html

# Add links to OBIGT.html
sed -i 's/reset()/<a href="..\/html\/thermo.html" style="color:\ red;">reset()<\/a>/g' OBIGT.html
# Note closing ) so we don't clash with add.OBIGT(
sed -i 's/OBIGT()/<a href="..\/html\/thermo.html" style="color:\ red;">OBIGT()<\/a>/g' OBIGT.html
sed -i 's/thermo()/<a href="..\/html\/thermo.html" style="color:\ green;">thermo()<\/a>/g' OBIGT.html
sed -i 's/add.OBIGT(/<a href="..\/html\/add.OBIGT.html" style="color:\ red;">add.OBIGT<\/a>(/g' OBIGT.html
sed -i 's/water(/<a href="..\/html\/water.html" style="color:\ green;">water<\/a>(/g' OBIGT.html
sed -i 's/demo(/<a href="..\/demo">demo<\/a>(/g' OBIGT.html

# Add links to multi-metal.html 20200716
sed -i 's/affinity()/<a href="..\/html\/affinity.html">affinity()<\/a>/g' multi-metal.html
sed -i 's/mash()/<a href="..\/html\/mix.html">mash()<\/a>/g' multi-metal.html
sed -i 's/diagram()/<a href="..\/html\/diagram.html">diagram()<\/a>/g' multi-metal.html
sed -i 's/stack_mosaic()/<a href="..\/html\/stack_mosaic.html">stack_mosaic()<\/a>/g' multi-metal.html
# Put a ">" in front so this doesn't match stack_mosaic 20230914
sed -i 's/>mosaic()/><a href="..\/html\/mosaic.html">mosaic()<\/a>/g' multi-metal.html
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
# Note closing ) so we don't clash with add.OBIGT(, mod.OBIGT(, or logK.to.OBIGT(
sed -i '375,$s/>OBIGT()/><a href="..\/html\/thermo.html" style="color:\ red;">OBIGT()<\/a>/g' custom_data.html
sed -i '375,$s/add.OBIGT(/<a href="..\/html\/add.OBIGT.html" style="color:\ red;">add.OBIGT<\/a>(/g' custom_data.html
sed -i '375,$s/mod.OBIGT(/<a href="..\/html\/add.OBIGT.html" style="color:\ red;">mod.OBIGT<\/a>(/g' custom_data.html
sed -i '375,$s/logK.to.OBIGT(/<a href="..\/html\/logK.to.OBIGT.html" style="color:\ red;">logK.to.OBIGT<\/a>(/g' custom_data.html
sed -i '375,$s/basis(/<a href="..\/html\/basis.html" style="color:\ red;">basis<\/a>(/g' custom_data.html
sed -i '375,$s/species(/<a href="..\/html\/species.html" style="color:\ red;">species<\/a>(/g' custom_data.html
sed -i '375,$s/E.units(/<a href="..\/html\/util.units.html" style="color:\ red;">E.units<\/a>(/g' custom_data.html
sed -i '375,$s/info(/<a href="..\/html\/info.html" style="color:\ green;">info<\/a>(/g' custom_data.html
sed -i '375,$s/subcrt(/<a href="..\/html\/subcrt.html" style="color:\ green;">subcrt<\/a>(/g' custom_data.html
sed -i '375,$s/affinity(/<a href="..\/html\/affinity.html" style="color:\ green;">affinity<\/a>(/g' custom_data.html
sed -i '375,$s/thermo.refs(/<a href="..\/html\/util.data.html" style="color:\ green;">thermo.refs<\/a>(/g' custom_data.html
sed -i '375,$s/thermo(/<a href="..\/html\/thermo.html" style="color:\ green;">thermo<\/a>(/g' custom_data.html
sed -i '375,$s/check.GHS(/<a href="..\/html\/util.data.html" style="color:\ green;">check.GHS<\/a>(/g' custom_data.html
sed -i '375,$s/check.EOS(/<a href="..\/html\/util.data.html" style="color:\ green;">check.EOS<\/a>(/g' custom_data.html

# Add links to FAQ.html 20230517
sed -i 's/<code>equilibrate()<\/code>/<code><a href="..\/html\/equilibrate.html" style="color:\ green;">equilibrate()<\/a><\/code>/g' FAQ.html
sed -i 's/<code>info()<\/code>/<code><a href="..\/html\/info.html" style="color:\ green;">info()<\/a><\/code>/g' FAQ.html
sed -i 's/<code>thermo.refs()<\/code>/<code><a href="..\/html\/util.data.html" style="color:\ green;">thermo.refs()<\/a><\/code>/g' FAQ.html
sed -i 's/<code>subcrt()<\/code>/<code><a href="..\/html\/subcrt.html" style="color:\ green;">subcrt()<\/a><\/code>/g' FAQ.html
sed -i 's/<code>check.GHS()<\/code>/<code><a href="..\/html\/util.data.html" style="color:\ green;">check.GHS()<\/a><\/code>/g' FAQ.html
sed -i 's/<code>affinity()<\/code>/<code><a href="..\/html\/affinity.html" style="color:\ green;">affinity()<\/a><\/code>/g' FAQ.html
sed -i 's/<code>diagram()<\/code>/<code><a href="..\/html\/diagram.html" style="color:\ green;">diagram()<\/a><\/code>/g' FAQ.html
sed -i 's/<code>mosaic()<\/code>/<code><a href="..\/html\/mosaic.html" style="color:\ green;">mosaic()<\/a><\/code>/g' FAQ.html
