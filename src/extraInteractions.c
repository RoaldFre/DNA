#include "extraInteractions.h"
#include "physics.h"

/* HARMONIC END-TO-END POTENTIAL */

static double harmonicEndToEndPotential(void *data)
{
	HarmonicEndToEndInt *hetei = (HarmonicEndToEndInt*) data;
	double R = endToEndDistance(hetei->s);
	return 0.5 * hetei->K * SQUARE(R - hetei->Rref);
}
static void harmonicEndToEndForce(void *data)
{
	HarmonicEndToEndInt *hetei = (HarmonicEndToEndInt*) data;
	Strand *s = hetei->s;
	double R = endToEndDistance(hetei->s);
	Vec3 Rdir = endToEndDirection(hetei->s);
	Vec3 F = scale(Rdir, hetei->K * (R - hetei->Rref));

	distributeForceOverMonomer(F,            s, 0);
	distributeForceOverMonomer(scale(F, -1), s, s->numMonomers - 1);
}

void registerHarmonicEndToEndInt(HarmonicEndToEndInt *conf)
{
	/* Register the end-to-end interaction */
	ExtraInteraction interaction;
	interaction.name = "endToEnd";
	interaction.symbol = "ete";
	interaction.data = conf;
	interaction.potential = &harmonicEndToEndPotential;
	interaction.addForces = &harmonicEndToEndForce;
	registerExtraInteraction(&interaction);
}
char *harmonicEndToEndIntHeader(HarmonicEndToEndInt *conf)
{
	return asprintfOrDie(
		"#\n"
		"# End to end interaction, with an umbrella sampling\n"
		"# potential on the first and last monomer of the strand.\n"
		"# This potential is given by\n"
		"#   V = K * (R - Rref),\n"
		"# with R the end-to-end distance between COM of first\n"
		"# and last monomer, and parameters:\n"
		"# <K> <Rref>\n"
		"%e\t%e\n"
		"#\n",
		conf->K, conf->Rref);
}

