let stepped = false;
document.addEventListener("DOMContentLoaded", async function () {
  const grabstep_normalization_rslt = await grabsteps(
    `/getSteps/normalization/process/`
  );
  const grabstep_merge_rslt = await grabsteps(`/getSteps/merge/process/`);

  if (
    grabstep_normalization_rslt.message != "notfound" &&
    grabstep_merge_rslt.message != "notfound"
  ) {
    stepped = true;
    document.getElementById("watchonly").classList.remove("hidden");
    document.getElementById("result").classList.remove("hidden");
    document.getElementById("method_select").disabled = true;
    show_result(grabstep_normalization_rslt.output_values, "n_adata_results");
    show_result(grabstep_merge_rslt.output_values, "m_adata_results");
    document.getElementById("nextbtn").classList.remove("hidden");
    banned_operations();
  }
});
function banned_operations() {
  const processBtn = document.getElementById("processbutton");
  processBtn.disabled = true;
}
// proccess button click
async function processbtn(event) {
  event.preventDefault();
  toggleLoading(true, "processbutton");
  await normalization();
  await merged();
  document.getElementById("result").classList.remove("hidden");
  // next page btn
  toggleLoading(false, "processbutton");
  document.getElementById("nextbtn").classList.remove("hidden");
}
async function normalization() {
  const csrftoken = getCookie("csrftoken");
  const form = document.getElementById("normalform");
  const formData = new FormData(form);
  const process_result = await fetchAPI("/api/normal", formData, csrftoken);
  // addata_results div
  show_result(process_result.data, "n_adata_results");
}
async function merged() {
  const csrftoken = getCookie("csrftoken");
  const process_result = await fetchAPI("/api/merge", 0, csrftoken);
  // addata_results div
  show_result(process_result.data, "m_adata_results");
}
function show_result(process_result, resultid) {
  let adata_results = document.getElementById(resultid);
  adata_results.innerHTML = "";
  let ul = document.createElement("ul");
  ul.classList.add("space-y-1", "text-gray-500", "list-disc", "list-inside");
  console.log(process_result["adata_results"]);
  process_result["adata_results"].forEach((item) => {
    let li = document.createElement("li");
    li.textContent = item;
    ul.appendChild(li);
  });
  adata_results.appendChild(ul);
}
function select(s_ul, selected) {
  document.getElementById(s_ul + "_select").textContent = selected;
  document.getElementById(s_ul + "_input").value = selected;
}
