import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html


lbt_ingest_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Scanpy Ingest Settings"),

        dbc.PopoverBody(
            [
                html.P(
                    [
                        "Details: ",
                        html.A("scanpy.readthedocs.io/en/"
                               "stable/api/scanpy.tl.ingest.html",
                               href="https://scanpy.readthedocs.io/en/"
                               "stable/api/scanpy.tl.ingest.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="lbt-Scanpy-Ingest-settings",
    target="lbt-Scanpy-Ingest-btn",
    trigger="click"
)

lbt_ingest_settings_keys = {

}

lbt_settings = [
    lbt_ingest_settings
]

lbt_settings_keys = {
    'lbt-Scanpy-Ingest': lbt_ingest_settings_keys,
}
