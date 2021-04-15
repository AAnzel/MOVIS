import { Axis as VgAxis, SignalRef, Text } from 'vega';
import { AxisInternal, AxisPart, AxisPropsWithCondition, ConditionalAxisProp } from '../../axis';
import { FieldDefBase } from '../../channeldef';
import { Split } from '../split';
export declare type AxisComponentProps = Omit<VgAxis, 'title' | ConditionalAxisProp> & Omit<AxisPropsWithCondition<SignalRef>, 'title'> & {
    title: Text | FieldDefBase<string>[] | SignalRef;
    labelExpr: string;
    disable: boolean;
};
export declare const AXIS_COMPONENT_PROPERTIES: ("title" | "values" | "orient" | "tickCount" | "aria" | "description" | "offset" | "titleAlign" | "titleAnchor" | "titleBaseline" | "titleColor" | "titleFont" | "titleFontSize" | "titleFontStyle" | "titleFontWeight" | "titleLimit" | "titleLineHeight" | "titleOpacity" | "titlePadding" | "labelAlign" | "labelBaseline" | "labelColor" | "labelFont" | "labelFontSize" | "labelFontStyle" | "labelFontWeight" | "labelLimit" | "labelOpacity" | "labelPadding" | "labelOffset" | "labelOverlap" | "labelSeparation" | "zindex" | "scale" | "domain" | "ticks" | "gridColor" | "gridDash" | "gridDashOffset" | "gridOpacity" | "gridWidth" | "tickColor" | "tickDash" | "tickDashOffset" | "tickOpacity" | "tickSize" | "tickWidth" | "gridScale" | "format" | "formatType" | "position" | "tickMinStep" | "encode" | "translate" | "minExtent" | "maxExtent" | "bandPosition" | "titleAngle" | "titleX" | "titleY" | "domainCap" | "domainDash" | "domainDashOffset" | "domainColor" | "domainOpacity" | "domainWidth" | "tickBand" | "tickCap" | "tickExtra" | "tickOffset" | "tickRound" | "grid" | "gridCap" | "labels" | "labelBound" | "labelFlush" | "labelFlushOffset" | "labelLineHeight" | "labelAngle" | "labelExpr" | "disable")[];
export declare class AxisComponent extends Split<AxisComponentProps> {
    readonly explicit: Partial<AxisComponentProps>;
    readonly implicit: Partial<AxisComponentProps>;
    mainExtracted: boolean;
    constructor(explicit?: Partial<AxisComponentProps>, implicit?: Partial<AxisComponentProps>, mainExtracted?: boolean);
    clone(): AxisComponent;
    hasAxisPart(part: AxisPart): boolean;
    hasOrientSignalRef(): boolean;
}
export interface AxisComponentIndex {
    x?: AxisComponent[];
    y?: AxisComponent[];
}
export interface AxisInternalIndex {
    x?: AxisInternal;
    y?: AxisInternal;
}
//# sourceMappingURL=component.d.ts.map