var __rest = (this && this.__rest) || function (s, e) {
    var t = {};
    for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p) && e.indexOf(p) < 0)
        t[p] = s[p];
    if (s != null && typeof Object.getOwnPropertySymbols === "function")
        for (var i = 0, p = Object.getOwnPropertySymbols(s); i < p.length; i++) {
            if (e.indexOf(p[i]) < 0 && Object.prototype.propertyIsEnumerable.call(s, p[i]))
                t[p[i]] = s[p[i]];
        }
    return t;
};
import { isArray } from 'vega';
import { isBinParams } from '../bin';
import { isConditionalDef, isFieldDef, isScaleFieldDef } from '../channeldef';
import { normalizeLogicalComposition } from '../logical';
import { SpecMapper } from '../spec/map';
import { isBin, isFilter, isLookup } from '../transform';
import { duplicate, entries, vals } from '../util';
export class SelectionCompatibilityNormalizer extends SpecMapper {
    map(spec, normParams) {
        var _a, _b;
        (_a = normParams.emptySelections) !== null && _a !== void 0 ? _a : (normParams.emptySelections = {});
        (_b = normParams.selectionPredicates) !== null && _b !== void 0 ? _b : (normParams.selectionPredicates = {});
        spec = normalizeTransforms(spec, normParams);
        return super.map(spec, normParams);
    }
    mapLayerOrUnit(spec, normParams) {
        spec = normalizeTransforms(spec, normParams);
        if (spec.encoding) {
            const encoding = {};
            for (const [channel, enc] of entries(spec.encoding)) {
                encoding[channel] = normalizeChannelDef(enc, normParams);
            }
            spec = Object.assign(Object.assign({}, spec), { encoding });
        }
        return super.mapLayerOrUnit(spec, normParams);
    }
    mapUnit(spec, normParams) {
        const _a = spec, { selection } = _a, rest = __rest(_a, ["selection"]);
        if (selection) {
            return Object.assign(Object.assign({}, rest), { params: entries(selection).map(([name, selDef]) => {
                    var _a;
                    const _b = selDef, { init: value, bind, empty } = _b, select = __rest(_b, ["init", "bind", "empty"]);
                    if (select.type === 'single') {
                        select.type = 'point';
                        select.toggle = false;
                    }
                    else if (select.type === 'multi') {
                        select.type = 'point';
                    }
                    // Propagate emptiness forwards and backwards
                    normParams.emptySelections[name] = empty !== 'none';
                    for (const pred of vals((_a = normParams.selectionPredicates[name]) !== null && _a !== void 0 ? _a : {})) {
                        pred.empty = empty !== 'none';
                    }
                    return { name, value, select, bind };
                }) });
        }
        return spec;
    }
}
function normalizeTransforms(spec, normParams) {
    const { transform: tx } = spec, rest = __rest(spec, ["transform"]);
    if (tx) {
        const transform = tx.map((t) => {
            if (isFilter(t)) {
                return { filter: normalizePredicate(t, normParams) };
            }
            else if (isBin(t) && isBinParams(t.bin)) {
                return Object.assign(Object.assign({}, t), { bin: normalizeBinExtent(t.bin) });
            }
            else if (isLookup(t)) {
                const _a = t.from, { selection: param } = _a, from = __rest(_a, ["selection"]);
                return param
                    ? Object.assign(Object.assign({}, t), { from: Object.assign({ param }, from) }) : t;
            }
            return t;
        });
        return Object.assign(Object.assign({}, rest), { transform });
    }
    return spec;
}
function normalizeChannelDef(obj, normParams) {
    var _a, _b;
    const enc = duplicate(obj);
    if (isFieldDef(enc) && isBinParams(enc.bin)) {
        enc.bin = normalizeBinExtent(enc.bin);
    }
    if (isScaleFieldDef(enc) && ((_b = (_a = enc.scale) === null || _a === void 0 ? void 0 : _a.domain) === null || _b === void 0 ? void 0 : _b.selection)) {
        const _c = enc.scale.domain, { selection: param } = _c, domain = __rest(_c, ["selection"]);
        enc.scale.domain = Object.assign(Object.assign({}, domain), (param ? { param } : {}));
    }
    if (isConditionalDef(enc)) {
        if (isArray(enc.condition)) {
            enc.condition = enc.condition.map((c) => {
                const { selection, param, test } = c, cond = __rest(c, ["selection", "param", "test"]);
                return param ? c : Object.assign(Object.assign({}, cond), { test: normalizePredicate(c, normParams) });
            });
        }
        else {
            const _d = normalizeChannelDef(enc.condition, normParams), { selection, param, test } = _d, cond = __rest(_d, ["selection", "param", "test"]);
            enc.condition = param
                ? enc.condition
                : Object.assign(Object.assign({}, cond), { test: normalizePredicate(enc.condition, normParams) });
        }
    }
    return enc;
}
function normalizeBinExtent(bin) {
    const ext = bin.extent;
    if (ext === null || ext === void 0 ? void 0 : ext.selection) {
        const { selection: param } = ext, rest = __rest(ext, ["selection"]);
        return Object.assign(Object.assign({}, bin), { extent: Object.assign(Object.assign({}, rest), { param }) });
    }
    return bin;
}
function normalizePredicate(op, normParams) {
    // Normalize old compositions of selection names (e.g., selection: {and: ["one", "two"]})
    const normalizeSelectionComposition = (o) => {
        return normalizeLogicalComposition(o, param => {
            var _a, _b;
            var _c;
            const empty = (_a = normParams.emptySelections[param]) !== null && _a !== void 0 ? _a : true;
            const pred = { param, empty };
            (_b = (_c = normParams.selectionPredicates)[param]) !== null && _b !== void 0 ? _b : (_c[param] = []);
            normParams.selectionPredicates[param].push(pred);
            return pred;
        });
    };
    return op.selection
        ? normalizeSelectionComposition(op.selection)
        : normalizeLogicalComposition(op.test || op.filter, o => o.selection ? normalizeSelectionComposition(o.selection) : o);
}
//# sourceMappingURL=selectioncompat.js.map