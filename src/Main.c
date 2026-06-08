#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef int16_t  i16;
typedef int32_t  i32;

#define READ_BE16(p) (((u16)(p)[0] << 8) | (p)[1])
#define READ_BE32(p) (((u32)(p)[0] << 24) | ((u32)(p)[1] << 16) | ((u32)(p)[2] << 8) | (p)[3])

// ====================== Strukturen ======================

typedef struct {
    u16 numberOfContours;
    i16 xMin, yMin, xMax, yMax;
} GlyphHeader;

typedef struct {
    u16 *endPtsOfContours;
    u16 instructionLength;
    u8 *flags;
    u16 numPoints;
} SimpleGlyph;

typedef struct {
    GlyphHeader header;
    bool isSimple;
    SimpleGlyph simple;
} GlyphData;

typedef struct {
    u8* data;
    size_t size;
    u16 numGlyphs;
    u16 indexToLocFormat;
    u32 glyfOffset;
    u8* locaData;

    // cmap Format 4
    u16 cmapSegCount;
    u16* endCode;
    u16* startCode;
    i16* idDelta;
    u16* idRangeOffset;
    u16* glyphIdArray;
} TTF_Font;

// ====================== Hilfsfunktionen ======================

static u32 find_table_offset(TTF_Font* font, const char* tag) {
    u16 numTables = READ_BE16(font->data + 4);
    u32 dir = 12;
    for (u16 i = 0; i < numTables; i++) {
        u32 off = dir + i * 16;
        if (memcmp(font->data + off, tag, 4) == 0) {
            return READ_BE32(font->data + off + 8);
        }
    }
    return 0;
}

int ttf_load(TTF_Font* font, const char* filename) {
    memset(font, 0, sizeof(*font));
    
    FILE* f = fopen(filename, "rb");
    if (!f) return 0;
    fseek(f, 0, SEEK_END);
    font->size = ftell(f);
    rewind(f);
    font->data = malloc(font->size);
    if (!font->data) { fclose(f); return 0; }
    fread(font->data, 1, font->size, f);
    fclose(f);

    // head
    u32 off = find_table_offset(font, "head");
    if (off) font->indexToLocFormat = READ_BE16(font->data + off + 50);

    // maxp
    off = find_table_offset(font, "maxp");
    if (off) font->numGlyphs = READ_BE16(font->data + off + 4);

    // loca & glyf
    off = find_table_offset(font, "loca");
    if (off) font->locaData = font->data + off;
    font->glyfOffset = find_table_offset(font, "glyf");

    // cmap Format 4
    off = find_table_offset(font, "cmap");
    if (off) {
        u16 numSub = READ_BE16(font->data + off + 2);
        for (u16 i = 0; i < numSub; i++) {
            u32 sub = off + 4 + i*8;
            u16 plat = READ_BE16(font->data + sub);
            u16 enc = READ_BE16(font->data + sub + 2);
            u32 cmapOff = READ_BE32(font->data + sub + 4);
            if ((plat == 3 && enc == 1) || plat == 0) {
                u32 fmtOff = off + cmapOff;
                if (READ_BE16(font->data + fmtOff) == 4) {
                    u16 segX2 = READ_BE16(font->data + fmtOff + 6);
                    font->cmapSegCount = segX2 / 2;
                    u32 eOff = fmtOff + 14;
                    font->endCode = (u16*)(font->data + eOff);
                    font->startCode = font->endCode + font->cmapSegCount + 1;
                    font->idDelta = (i16*)(font->startCode + font->cmapSegCount);
                    font->idRangeOffset = (u16*)(font->idDelta + font->cmapSegCount);
                    font->glyphIdArray = font->idRangeOffset + font->cmapSegCount;
                    break;
                }
            }
        }
    }
    return 1;
}

u16 ttf_get_glyph_index(TTF_Font* font, u32 codepoint) {
    if (!font->endCode) return 0;
    for (int i = 0; i < font->cmapSegCount; i++) {
        u16 end = READ_BE16((u8*)&font->endCode[i]);
        if (codepoint > end) continue;
        u16 start = READ_BE16((u8*)&font->startCode[i]);
        if (codepoint < start) break;

        if (font->idRangeOffset[i] == 0) {
            return (u16)(codepoint + READ_BE16((u8*)&font->idDelta[i]));
        } else {
            //u16 off = READ_BE16((u8*)&font->idRangeOffset[i]) / 2 + (codepoint - start);
            //u16 gid = READ_BE16((u8*)&font->glyphIdArray[off]);

            u16 rangeOffset = READ_BE16((u8*)&font->idRangeOffset[i]) / 2 + (codepoint - start);
            u8* roPtr =(u8*)&font->idRangeOffset[i];
            u8* glyphPtr = roPtr + rangeOffset + 2 * (codepoint - start);
            u16 gid = READ_BE16(glyphPtr);
            return gid + READ_BE16((u8*)&font->idDelta[i]);
        }
    }
    return 0;
}

static u32 get_glyph_offset(TTF_Font* font, u16 glyphIndex) {
    if (glyphIndex >= font->numGlyphs) return 0;
    if (font->indexToLocFormat == 0) {
        return READ_BE16(font->locaData + glyphIndex * 2) * 2;
    } else {
        return READ_BE32(font->locaData + glyphIndex * 4);
    }
}

// ====================== Vollständiges Punkt-Parsing ======================

typedef struct {
    i16 x, y;
    bool onCurve;
} Point;

void parse_simple_glyph_points(u8* glyphData, GlyphData* g, Point** ppoints) {
    u32 ptr = 10; // nach Header

    g->simple.endPtsOfContours = (u16*)(glyphData + ptr);
    ptr += g->header.numberOfContours * sizeof(u16);   // wichtig: header verwenden!

    g->simple.instructionLength = READ_BE16(glyphData + ptr);
    ptr += 2 + g->simple.instructionLength;

    u16 lastEnd = READ_BE16((u8*)&g->simple.endPtsOfContours[g->header.numberOfContours - 1]);
    g->simple.numPoints = lastEnd + 1;

    *ppoints = malloc(g->simple.numPoints * sizeof(Point));
    Point* points = *ppoints;

    g->simple.flags = malloc(g->simple.numPoints);

    // Flags + REPEAT
    u8* flagStart = glyphData + ptr;
    u16 flagIdx = 0, ptIdx = 0;

    while (ptIdx < g->simple.numPoints) {
        u8 flag = flagStart[flagIdx++];
        u16 repeat = 1;
        if (flag & 0x08) repeat += flagStart[flagIdx++];

        for (u16 r = 0; r < repeat && ptIdx < g->simple.numPoints; r++) {
            g->simple.flags[ptIdx] = flag;
            points[ptIdx].onCurve = (flag & 1) != 0;
            ptIdx++;
        }
    }
    ptr += flagIdx;

    // X-Koordinaten
    i16 x = 0;
    for (u16 i = 0; i < g->simple.numPoints; i++) {
        u8 flag = g->simple.flags[i];
        if (flag & 0x02) {                         // X_SHORT
            if (flag & 0x10) x += glyphData[ptr++];
            else             x -= glyphData[ptr++];
        } else if (!(flag & 0x10)) {               // long delta
            x += READ_BE16(glyphData + ptr);
            ptr += 2;
        }
        points[i].x = x;
    }

    // Y-Koordinaten
    i16 y = 0;
    for (u16 i = 0; i < g->simple.numPoints; i++) {
        u8 flag = g->simple.flags[i];
        if (flag & 0x04) {                         // Y_SHORT
            if (flag & 0x20) y += glyphData[ptr++];
            else             y -= glyphData[ptr++];
        } else if (!(flag & 0x20)) {
            y += READ_BE16(glyphData + ptr);
            ptr += 2;
        }
        points[i].y = y;
    }
}

// ====================== Öffentliche Funktionen ======================

GlyphData ttf_load_glyph_full(TTF_Font* font, u16 glyphIndex, Point** out_points) {
    GlyphData g = {0};
    u32 offset0 = get_glyph_offset(font, glyphIndex);
    u32 offset1 = get_glyph_offset(font, glyphIndex + 1);

    //if (offset0 == offset1) return g;
    if (!offset0) return g;

    u8* glyphData = font->data + font->glyfOffset + offset0;

    g.header.numberOfContours = READ_BE16(glyphData);
    g.header.xMin = READ_BE16(glyphData + 2);
    g.header.yMin = READ_BE16(glyphData + 4);
    g.header.xMax = READ_BE16(glyphData + 6);
    g.header.yMax = READ_BE16(glyphData + 8);

    if (g.header.numberOfContours > 0) {
        g.isSimple = true;
        SimpleGlyph* s = &g.simple;

        Point* points = NULL;
        parse_simple_glyph_points(glyphData, &g, &points);
        *out_points = points;
    }
    return g;
}

void ttf_free_glyph(GlyphData* g, Point* points) {
    if (g->isSimple) {
        free(g->simple.flags);
        free(points);
    }
}

void ttf_free(TTF_Font* font) {
    free(font->data);
}

// ====================== main ======================
#define FONT_PATH "./assets/JetBrainsMono-Bold.ttf"

int main(int argc, char** argv) {
    TTF_Font font = {0};
    if (!ttf_load(&font, FONT_PATH)) {
        printf("Konnte Font nicht laden!\n");
        return 1;
    }

    printf("Font geladen: %d Glyphen\n", font.numGlyphs);

    u16 gid = ttf_get_glyph_index(&font, 'A');
    Point* points = NULL;
    GlyphData g = ttf_load_glyph_full(&font, gid, &points);

    if (g.isSimple && points) {
        printf("Glyph 'A': %d Konturen, %d Punkte\n", 
               g.header.numberOfContours, g.simple.numPoints);

        u16 start = 0;
        for (u16 c = 0; c < g.header.numberOfContours; c++) {
            u16 end = READ_BE16((u8*)&g.simple.endPtsOfContours[c]);
            printf("Kontur %d (Punkte %d-%d):\n", c, start, end);

            for (u16 i = start; i <= end; i++) {
                Point p = points[i];
                printf("  P%3d: (%5d, %5d) %s\n", i, p.x, p.y, p.onCurve ? "ON " : "OFF");
            }
            start = end + 1;
        }

        ttf_free_glyph(&g, points);
    }

    ttf_free(&font);
    return 0;
}