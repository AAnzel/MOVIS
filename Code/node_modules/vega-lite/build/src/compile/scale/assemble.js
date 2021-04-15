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
import { isObject } from 'vega-util';
import { isXorY } from '../../channel';
import { keys } from '../../util';
import { isDataRefDomain, isVgRangeStep } from '../../vega.schema';
import { isConcatModel, isLayerModel } from '../model';
import { assembleSelectionScaleDomain } from '../selection/assemble';
import { assembleDomain } from './domain';
export function assembleScales(model) {
    if (isLayerModel(model) || isConcatModel(model)) {
        // For concat and layer, include scales of children too
        return model.children.reduce((scales, child) => {
            return scales.concat(assembleScales(child));
        }, assembleScalesForModel(model));
    }
    else {
        // For facet, child scales would not be included in the parent's scope.
        // For unit, there is no child.
        return assembleScalesForModel(model);
    }
}
export function assembleScalesForModel(model) {
    return keys(model.component.scales).reduce((scales, channel) => {
        const scaleComponent = model.component.scales[channel];
        if (scaleComponent.merged) {
            // Skipped merged scales
            return scales;
        }
        const scale = scaleComponent.combine();
        const { name, type, selectionExtent, domains: _d, range: _r, reverse } = scale, otherScaleProps = __rest(scale, ["name", "type", "selectionExtent", "domains", "range", "reverse"]);
        const range = assembleScaleRange(scale.range, name, channel, model);
        const domain = assembleDomain(model, channel);
        const domainRaw = selectionExtent
            ? assembleSelectionScaleDomain(model, selectionExtent, scaleComponent, domain)
            : null;
        scales.push(Object.assign(Object.assign(Object.assign(Object.assign(Object.assign({ name,
            type }, (domain ? { domain } : {})), (domainRaw ? { domainRaw } : {})), { range }), (reverse !== undefined ? { reverse: reverse } : {})), otherScaleProps));
        return scales;
    }, []);
}
export function assembleScaleRange(scaleRange, scaleName, channel, model) {
    // add signals to x/y range
    if (isXorY(channel)) {
        if (isVgRangeStep(scaleRange)) {
            // For width/height step, use a signal created in layout assemble instead of a constant step.
            return {
                step: { signal: `${scaleName}_step` }
            };
        }
    }
    else if (isObject(scaleRange) && isDataRefDomain(scaleRange)) {
        return Object.assign(Object.assign({}, scaleRange), { data: model.lookupDataSource(scaleRange.data) });
    }
    return scaleRange;
}
//# sourceMappingURL=assemble.js.map