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
import { vgField } from '../../channeldef';
import { getTimeUnitParts, normalizeTimeUnit } from '../../timeunit';
import { duplicate, entries, hash, isEmpty, replacePathInField, vals } from '../../util';
import { DataFlowNode } from './dataflow';
export class TimeUnitNode extends DataFlowNode {
    constructor(parent, formula) {
        super(parent);
        this.formula = formula;
    }
    clone() {
        return new TimeUnitNode(null, duplicate(this.formula));
    }
    static makeFromEncoding(parent, model) {
        const formula = model.reduceFieldDef((timeUnitComponent, fieldDef) => {
            const { field, timeUnit } = fieldDef;
            if (timeUnit) {
                const as = vgField(fieldDef, { forAs: true });
                timeUnitComponent[hash({
                    as,
                    field,
                    timeUnit
                })] = {
                    as,
                    field,
                    timeUnit
                };
            }
            return timeUnitComponent;
        }, {});
        if (isEmpty(formula)) {
            return null;
        }
        return new TimeUnitNode(parent, formula);
    }
    static makeFromTransform(parent, t) {
        const _a = Object.assign({}, t), { timeUnit } = _a, other = __rest(_a, ["timeUnit"]);
        const normalizedTimeUnit = normalizeTimeUnit(timeUnit);
        const component = Object.assign(Object.assign({}, other), { timeUnit: normalizedTimeUnit });
        return new TimeUnitNode(parent, {
            [hash(component)]: component
        });
    }
    /**
     * Merge together TimeUnitNodes assigning the children of `other` to `this`
     * and removing `other`.
     */
    merge(other) {
        this.formula = Object.assign({}, this.formula);
        // if the same hash happen twice, merge
        for (const key in other.formula) {
            if (!this.formula[key]) {
                // copy if it's not a duplicate
                this.formula[key] = other.formula[key];
            }
        }
        for (const child of other.children) {
            other.removeChild(child);
            child.parent = this;
        }
        other.remove();
    }
    /**
     * Remove time units coming from the other node.
     */
    removeFormulas(fields) {
        const newFormula = {};
        for (const [key, timeUnit] of entries(this.formula)) {
            if (!fields.has(timeUnit.as)) {
                newFormula[key] = timeUnit;
            }
        }
        this.formula = newFormula;
    }
    producedFields() {
        return new Set(vals(this.formula).map(f => f.as));
    }
    dependentFields() {
        return new Set(vals(this.formula).map(f => f.field));
    }
    hash() {
        return `TimeUnit ${hash(this.formula)}`;
    }
    assemble() {
        const transforms = [];
        for (const f of vals(this.formula)) {
            const { field, as, timeUnit } = f;
            const _a = normalizeTimeUnit(timeUnit), { unit, utc } = _a, params = __rest(_a, ["unit", "utc"]);
            transforms.push(Object.assign(Object.assign(Object.assign(Object.assign({ field: replacePathInField(field), type: 'timeunit' }, (unit ? { units: getTimeUnitParts(unit) } : {})), (utc ? { timezone: 'utc' } : {})), params), { as: [as, `${as}_end`] }));
        }
        return transforms;
    }
}
//# sourceMappingURL=timeunit.js.map