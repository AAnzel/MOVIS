import { isArray, isString } from 'vega';
import { isSelectionParameter } from '../selection';
import { isUnitSpec } from '../spec';
import { SpecMapper } from '../spec/map';
export class TopLevelSelectionsNormalizer extends SpecMapper {
    map(spec, normParams) {
        var _a;
        const selections = (_a = normParams.selections) !== null && _a !== void 0 ? _a : [];
        if (spec.params && !isUnitSpec(spec)) {
            const params = [];
            for (const param of spec.params) {
                if (isSelectionParameter(param)) {
                    selections.push(param);
                }
                else {
                    params.push(param);
                }
            }
            spec.params = params;
        }
        normParams.selections = selections;
        return super.map(spec, addSpecNameToParams(spec, normParams));
    }
    mapUnit(spec, normParams) {
        var _a;
        const selections = normParams.selections;
        if (!selections || !selections.length)
            return spec;
        const path = ((_a = normParams.path) !== null && _a !== void 0 ? _a : []).concat(spec.name);
        const params = [];
        for (const selection of selections) {
            // By default, apply selections to all unit views.
            if (!selection.views || !selection.views.length) {
                params.push(selection);
            }
            else {
                for (const view of selection.views) {
                    // view is either a specific unit name, or a partial path through the spec tree.
                    if ((isString(view) && (view === spec.name || path.indexOf(view) >= 0)) ||
                        (isArray(view) &&
                            view.map(v => path.indexOf(v)).every((v, i, arr) => v !== -1 && (i === 0 || v > arr[i - 1])))) {
                        params.push(selection);
                    }
                }
            }
        }
        if (params.length)
            spec.params = params;
        return spec;
    }
}
for (const method of ['mapFacet', 'mapRepeat', 'mapHConcat', 'mapVConcat', 'mapLayer']) {
    const proto = TopLevelSelectionsNormalizer.prototype[method];
    TopLevelSelectionsNormalizer.prototype[method] = function (spec, params) {
        return proto.call(this, spec, addSpecNameToParams(spec, params));
    };
}
function addSpecNameToParams(spec, params) {
    var _a;
    return spec.name
        ? Object.assign(Object.assign({}, params), { path: ((_a = params.path) !== null && _a !== void 0 ? _a : []).concat(spec.name) }) : params;
}
//# sourceMappingURL=toplevelselection.js.map