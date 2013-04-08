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
	ExtraInteraction interaction;
	interaction.name = "harmonicEndToEnd";
	interaction.symbol = "hete";
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
		"#   V = K/2 * (R - Rref)^2,\n"
		"# with R the end-to-end distance between COM of first\n"
		"# and last monomer, and parameters:\n"
		"# K = %le\n"
		"# Rref = %le\n"
		"#\n",
		conf->K, conf->Rref);
}



/* CONSTANT END-TO-END-FORCE */

static double endToEndForcePotential(void *data)
{
	EndToEndForceInt *etefi = (EndToEndForceInt*) data;
	Vec3 R = endToEndVector(etefi->s);
	return -dot(R, etefi->F);
}

static void endToEndForceForce(void *data)
{
	EndToEndForceInt *etefi = (EndToEndForceInt*) data;
	Vec3 F = etefi->F;
	Strand *s = etefi->s;
	distributeForceOverMonomer(scale(F, -1), s, 0);
	distributeForceOverMonomer(F,            s, s->numMonomers - 1);
}

void registerEndToEndForceInt(EndToEndForceInt *conf)
{
	ExtraInteraction interaction;
	interaction.name = "constantEndToEnd";
	interaction.symbol = "cete";
	interaction.data = conf;
	interaction.potential = &endToEndForcePotential;
	interaction.addForces = &endToEndForceForce;
	registerExtraInteraction(&interaction);
}

char *endToEndForceIntHeader(EndToEndForceInt *conf)
{
	return asprintfOrDie(
		"#\n"
		"# Constant end-to-end force: -F on first monomer, +F on last monomer.\n"
		"# F = %le %le %le\n"
		"#\n",
		conf->F.x, conf->F.y, conf->F.z);
}

