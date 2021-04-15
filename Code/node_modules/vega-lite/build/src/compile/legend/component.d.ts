import { Legend as VgLegend } from 'vega';
import { NonPositionScaleChannel } from '../../channel';
import { LegendInternal } from '../../legend';
import { Split } from '../split';
export declare type LegendComponentProps = VgLegend & {
    labelExpr?: string;
    selections?: string[];
    disable?: boolean;
};
export declare const LEGEND_COMPONENT_PROPERTIES: ("labelExpr" | "disable" | keyof VgLegend | "selections")[];
export declare class LegendComponent extends Split<LegendComponentProps> {
}
export declare type LegendComponentIndex = Partial<Record<NonPositionScaleChannel, LegendComponent>>;
export declare type LegendInternalIndex = Partial<Record<NonPositionScaleChannel, LegendInternal>>;
//# sourceMappingURL=component.d.ts.map