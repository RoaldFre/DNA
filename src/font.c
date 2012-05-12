#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdarg.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_GLYPH_H

#include "font.h"
static char *strdup(const char *string)
{
	char *new;
	new = malloc(strlen(string) + 1);
	strcpy(new, string);
	return new;
}

static FT_Library fontlib = NULL;

static void fontlib_destroy(void)
{
	if (fontlib)
		FT_Done_FreeType(fontlib);
}

Font *font_load(const char *filename)
{
	Font *font;

	if (fontlib == NULL)
	{
		if (FT_Init_FreeType(&fontlib) != 0)
		{
			printf("Error initializing FreeType 2\n");
			return NULL;
		}
		atexit(fontlib_destroy);
	}
	if ((font = malloc(sizeof(Font))) == NULL)
	{
		printf("Out of memory\n");
		return NULL;
	}
	if (FT_New_Face(fontlib, filename, 0, &font->face) != 0)
	{
		printf("Error initializing font %s\n", filename);
		free(font);
		return NULL;
	}
	printf("Loaded font with %ld glyphs\n", font->face->num_glyphs);

	return font;
}

void font_destroy(Font *font)
{
	FT_Done_Face(font->face);
	free(font);
}

static const unsigned int utf8_tail_length[256] = {
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x0F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x1F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x2F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x3F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x4F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x5F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x6F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x7F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x8F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0x9F */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0xAF */
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* 0xBF */
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 0xCF */
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, /* 0xDF */
2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, /* 0xEF */
3,3,3,3,3,3,3,3,4,4,4,4,5,5,0,0, /* 0xFF */
};

static uint32_t utf8_next(const uint8_t **s)
{
	uint32_t c;
	int i, tail_len;

	/* Start with the first byte */
	c = (*s)[0];
	(*s)++;

	/* Return early if it's plain ASCII */
	if (c < 0x80)
		return c;

	/* Store the first 1-6 bits */
	tail_len = utf8_tail_length[c & 0xff];
	c &= (0x3f >> tail_len);

	/* Main decoding loop */
	for (i = 0; i < tail_len; i++)
	{
		if ((*s)[i] == '\0' || ((*s)[i] & 0xc0) != 0x80)
			break; /* End of string or invalid continuation byte */

		c = (c << 6) + ((*s)[i] & 0x3f);
	}

	*s += i;
	if (i != tail_len)
		return 0xFFFD; /* Replacement character � */

	/* TODO: Check for overlong encodings, surrogate pairs etc */

	return c;
}

static void glyphstring_create(FT_Face face, Text *text, FT_Glyph *glyph_string,
		FT_Vector *pos)
{
	const uint8_t *string = text->string;
	FT_Bool has_kerning;
	FT_UInt glyph_index, previous;
	FT_Vector pen, delta;
	uint32_t charcode;
	int i;

	has_kerning = FT_HAS_KERNING(face);
	previous = 0;
	i = 0;
	pen.x = pen.y = 0;
	while (string[0] != '\0')
	{
		charcode = utf8_next(&string);
		glyph_index = FT_Get_Char_Index(face, charcode);
		if (has_kerning && previous && glyph_index)
			FT_Get_Kerning(face, previous, glyph_index, FT_KERNING_DEFAULT,
					&delta);
		else
			delta.x = 0;

		if (glyph_index == 0)
			printf("Glyph for character U+%X missing\n", charcode);

		if (FT_Load_Glyph(face, glyph_index, FT_LOAD_RENDER) != 0)
		{
			printf("Error loading glyph for character U+%X\n", charcode);
			continue;
		}
		if (FT_Get_Glyph(face->glyph, &glyph_string[i]) != 0)
		{
			printf("Error copying glyph for character U+%X\n", charcode);
			continue;
		}

		pen.x += delta.x;

		pos[i] = pen;

		pen.x += face->glyph->advance.x;
		pen.y += face->glyph->advance.y;

		previous = glyph_index;
		i++;
	}

	text->num_glyphs = i;
}

static void compute_glyphstring_bbox(FT_Glyph *string, FT_Vector *pos,
		int len, FT_BBox *bbox)
{
	FT_BBox glyph_bbox;
	int i;

	bbox->xMin = bbox->yMin = 32000;
	bbox->xMax = bbox->yMax = -32000;

	for (i = 0; i < len; i++)
	{
		FT_Glyph_Get_CBox(string[i], FT_GLYPH_BBOX_GRIDFIT, &glyph_bbox);

		glyph_bbox.xMin += pos[i].x;
		glyph_bbox.xMax += pos[i].x;
		glyph_bbox.yMin += pos[i].y;
		glyph_bbox.yMax += pos[i].y;

		if (glyph_bbox.xMin < bbox->xMin)
			bbox->xMin = glyph_bbox.xMin;

		if (glyph_bbox.xMax > bbox->xMax)
			bbox->xMax = glyph_bbox.xMax;

		if (glyph_bbox.yMin < bbox->yMin)
			bbox->yMin = glyph_bbox.yMin;

		if (glyph_bbox.yMax > bbox->yMax)
			bbox->yMax = glyph_bbox.yMax;
	}

	if (bbox->xMin > bbox->xMax || bbox->yMin > bbox->yMax)
	{
		bbox->xMin = 0;
		bbox->xMax = 0;
		bbox->yMin = 0;
		bbox->yMax = 0;
	}
}

static void blit_glyph(FT_BitmapGlyph glyph, GLubyte *image, int x, int y,
		int width, int height)
{
	FT_Bitmap *bitmap = &glyph->bitmap;
	int startx, starty, imagex, imagey, ix, iy;

	/* The bitmap doesn't contain the complete glyph, only the nonzero bytes */
	startx = x + glyph->left;
	starty = y + glyph->top - bitmap->rows;

	for (iy = 0; iy < bitmap->rows; iy++)
	{
		for (ix = 0; ix < bitmap->width; ix++)
		{
			/* The glyph is stored upside-down */
			GLubyte b = bitmap->buffer[(bitmap->rows - 1 - iy)*bitmap->pitch
					+ ix];
			imagex = startx + ix;
			imagey = starty + iy;

			/* This won't happen but let's make sure */
			if (imagex < 0 || imagex > width || imagey < 0 || imagey > height)
				continue;

			if (image[imagey*width + imagex] == 0)
				image[imagey*width + imagex] =  b;
		}
	}
}

Text *text_create(Font *font, const char *char_string, int size)
{
	const uint8_t *string = (const uint8_t *) char_string;
	Text *text;
	FT_Face face = font->face;
	FT_Glyph *glyph_string;
	FT_Vector *pos;
	FT_BBox bbox;
	FT_Long width, height;
	size_t len;
	float x, y;
	int i;

	if (FT_Set_Char_Size(face, 0, size*64, 0, 0) != 0)
	{
		printf("Error setting character size\n");
		return NULL;
	}

	text = malloc(sizeof(Text));
	if (text == NULL)
	{
		printf("Out of memory\n");
		return NULL;
	}
	text->size = size;
	text->vao = text->vbo = text->texture = 0;
	text->texture_image = NULL;
	text->string = (uint8_t *) strdup((const char *) string);
	len = strlen((const char *) string); /* Bytecount */
	/* We allocate more space than necessary for the glyph string, better safe
	 * than sorry. */
	glyph_string = calloc(sizeof(FT_Glyph), len);
	pos = calloc(sizeof(FT_Vector), len);
	if (text->string == NULL || glyph_string == NULL || pos == NULL)
	{
		printf("Out of memory\n");
		free(glyph_string);
		free(text);
		return NULL;
	}

	/* The UTF-8 bytestring is converted to an array of glyphs */
	glyphstring_create(face, text, glyph_string, pos);
	/* We determine how big the text will be. This information is used
	 * to compute the size of the texture we'll store the text in */
	compute_glyphstring_bbox(glyph_string, pos, text->num_glyphs, &bbox);

	width  = bbox.xMax - bbox.xMin;
	height = bbox.yMax - bbox.yMin;
	text->width = (width/64 + 0x3) & ~0x3; /* Align to 4 bytes */
	text->height = height/64;
	text->texture_image = calloc(sizeof(GLubyte),
			text->width * text->height);
	if (text->texture_image == NULL)
	{
		printf("Out of memory\n");
		for (i = 0; i < text->num_glyphs; i++)
			FT_Done_Glyph(glyph_string[i]);
		free(glyph_string);
		free(text);
		return NULL;
	}
	/* Now we can render the text to a texture */
	for (i = 0; i < text->num_glyphs; i++)
	{
		FT_Glyph glyph;
		FT_BitmapGlyph bitmap_glyph;
		FT_Vector pen;

		glyph = glyph_string[i];
		pen.x = pos[i].x;
		pen.y = pos[i].y;

		/* Render the new glyph and destroy the old one */
		if (FT_Glyph_To_Bitmap(&glyph, FT_LOAD_TARGET_NORMAL, &pen, 1) != 0)
		{
			printf("Error rendering glyph to bitmap for character '%c'\n",
					string[i]);
			FT_Done_Glyph(glyph_string[i]);
			continue;
		}

		bitmap_glyph = (FT_BitmapGlyph) glyph;
		blit_glyph(bitmap_glyph, text->texture_image, (pen.x - bbox.xMin)/64,
				(pen.y - bbox.yMin)/64, text->width, text->height);
		FT_Done_Glyph(glyph);
	}
	free(glyph_string);
	free(pos);

	x = bbox.xMin/64;
	y = bbox.yMin/64;

	/* Why do we add text->width and not bbox.xMax? Because OpenGL wants
	 * textures aligned at 4 bytes. This fixes a very subtle bug with 
	 * a very noticeable effect */
	text->vertex[0].x = x;
	text->vertex[0].y = y;
	text->vertex[0].u = 0;
	text->vertex[0].v = 0;

	text->vertex[1].x = x + text->width;
	text->vertex[1].y = y;
	text->vertex[1].u = 1;
	text->vertex[1].v = 0;

	text->vertex[2].x = x + text->width;
	text->vertex[2].y = y + text->height;
	text->vertex[2].u = 1;
	text->vertex[2].v = 1;

	text->vertex[3].x = x;
	text->vertex[3].y = y + text->height;
	text->vertex[3].u = 0;
	text->vertex[3].v = 1;

	for (i = 0; i < 4; i++)
	{
		text->vertex[i].r = 0;
		text->vertex[i].g = 0;
		text->vertex[i].b = 1;
		text->vertex[i].a = 1;
	}

	return text;
}

void text_upload_to_gpu(Text *text)
{
	glGenTextures(1, &text->texture);
	glBindTexture(GL_TEXTURE_2D, text->texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, text->width, text->height, 0,
			GL_ALPHA, GL_UNSIGNED_BYTE, text->texture_image);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	free(text->texture_image);
}

void text_render(Text *text)
{
	glBindTexture(GL_TEXTURE_2D, text->texture);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex2f(text->vertex[0].x, text->vertex[0].y);
	glTexCoord2f(1.0, 0.0);
	glVertex2f(text->vertex[1].x, text->vertex[1].y);
	glTexCoord2f(1.0, 1.0);
	glVertex2f(text->vertex[2].x, text->vertex[2].y);
	glTexCoord2f(0.0, 1.0);
	glVertex2f(text->vertex[3].x, text->vertex[3].y);
	glEnd();
}

void text_destroy(Text *text)
{
	glDeleteTextures(1, &text->texture);
	free(text->string);
	free(text);
}

void text_create_and_render(Font *font, int size, const char *string)
{
	Text *text;

	if (string == NULL)
		return;

	text = text_create(font, string, size);
	if (text ==	NULL)
		return;

	text_upload_to_gpu(text);
	text_render(text);
	text_destroy(text);
}
