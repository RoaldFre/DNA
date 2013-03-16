#include "extraInteractions.h"
#include "physics.h"

static double endToEndPotential(void *data)
{
	EndToEndInteraction *etei = (EndToEndInteraction*) data;
	double R = endToEndDistance(etei->s);
	return 0.5 * etei->K * SQUARE(R - etei->Rref);
}
static void endToEndForce(void *data)
{
	EndToEndInteraction *etei = (EndToEndInteraction*) data;
	Strand *s = etei->s;
	double R = endToEndDistance(etei->s);
	Vec3 Rdir = endToEndDirection(etei->s);
	Vec3 F = scale(Rdir, etei->K * (R - etei->Rref));

	distributeForceOverMonomer(F,            s, 0);
	distributeForceOverMonomer(scale(F, -1), s, s->numMonomers - 1);
}

void registerEndToEndInteraction(EndToEndInteraction *conf)
{
	/* Register the end-to-end interaction */
	ExtraInteraction interaction;
	interaction.name = "endToEnd";
	interaction.symbol = "ete";
	interaction.data = conf;
	interaction.potential = &endToEndPotential;
	interaction.addForces = &endToEndForce;
	registerExtraInteraction(&interaction);
}
char *endToEndInteractionHeader(EndToEndInteraction *conf)
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
