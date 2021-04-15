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
import { duplicate, hash } from '../../util';
import { DataFlowNode } from './dataflow';
/**
 * A class for regression transform nodes
 */
export class RegressionTransformNode extends DataFlowNode {
    constructor(parent, transform) {
        var _a, _b, _c;
        super(parent);
        this.transform = transform;
        this.transform = duplicate(transform); // duplicate to prevent side effects
        const specifiedAs = (_a = this.transform.as) !== null && _a !== void 0 ? _a : [undefined, undefined];
        this.transform.as = [(_b = specifiedAs[0]) !== null && _b !== void 0 ? _b : transform.on, (_c = specifiedAs[1]) !== null && _c !== void 0 ? _c : transform.regression];
    }
    clone() {
        return new RegressionTransformNode(null, duplicate(this.transform));
    }
    dependentFields() {
        var _a;
        return new Set([this.transform.regression, this.transform.on, ...((_a = this.transform.groupby) !== null && _a !== void 0 ? _a : [])]);
    }
    producedFields() {
        return new Set(this.transform.as);
    }
    hash() {
        return `RegressionTransform ${hash(this.transform)}`;
    }
    assemble() {
        const _a = this.transform, { regression, on } = _a, rest = __rest(_a, ["regression", "on"]);
        const result = Object.assign({ type: 'regression', x: on, y: regression }, rest);
        return result;
    }
}
//# sourceMappingURL=regression.js.map